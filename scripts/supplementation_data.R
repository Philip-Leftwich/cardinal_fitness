source(here::here("scripts", "packages.R"))
source(here::here("scripts", "functions", "supplementation_data_functions.R"))
source(here::here("scripts", "functions", "blood_feeding_functions.R"))

# Load data ====

path <- here::here("data", "XASupplementationData.xlsx")

data <- path |>
  excel_sheets() |>
  set_names() |>
  map(\(s) read_excel(path = path, sheet = s)) |>
  list_rbind(names_to = "experiment") |>
  filter(experiment != "Additional assays")

# Shape data ====
data_long <- data |>
  pivot_longer(
    cols = `WT (0mM)`:`Het (24mM)`,
    names_to = "treatment",
    values_to = "event"
  ) |>
  drop_na(event)

data_long


# Format data for interaction models
# Note: the WT NaOH-control column is headed "WT(6mM C-)" (no space) in both
# sheets. It is normalised to "WT (6mM C-)" here; previously the split on
# " (" failed for this column and the WT vehicle controls were silently
# dropped by the dosage filter (n = 211 instead of 244).
data_long2 <- data_long |>
  mutate(treatment = str_replace(treatment, fixed("WT("), "WT (")) |>
  separate_wider_delim(
    treatment,
    delim = " (",
    names = c("genotype", "dosage"),
    too_few = "align_start"
  ) |>
  mutate(genotype = str_trim(genotype)) |>
  mutate(dosage = str_remove(dosage, fixed(")"))) |>
  filter(dosage %in% c("0mM", "6mM", "6mM C-")) |>
  mutate(experiment = str_remove(experiment, "XA ")) |>
  mutate(
    # "C-" is the 0.1 N NaOH vehicle control for the 6 mM XA stock
    dosage_mM = factor(dosage, levels = c("0mM", "6mM C-", "6mM")) |>
      fct_recode("C-" = "6mM C-")
  )

# Biological plausibility checks ----
# Hard checks: survival times must be positive, event indicator must be
# binary, and genotype/experiment must belong to expected categories;
# validated before any modelling (produces Figure 4)
data_long2 <- data_long2 |>
  verify(Hours > 0) |> # survival time must be positive
  verify(event %in% c(0, 1)) |> # event must be binary (0/1)
  verify(!is.na(dosage_mM)) |> # every row maps to 0mM, C- or 6mM
  assert(in_set("WT", "Het"), genotype) |> # only two expected genotypes
  assert(in_set("sugar supplementation", "blood supplementation"), experiment) # only two supplementation routes

# Every arm x genotype x treatment group is present (catches dropped columns)
stopifnot(nrow(distinct(data_long2, experiment, genotype, dosage_mM)) == 12)

# Observed deaths per arm x genotype x treatment
xa_observed <- data_long2 |>
  summarise(
    n = n(),
    deaths = sum(event),
    .by = c(experiment, genotype, dosage_mM)
  ) |>
  arrange(experiment, genotype, dosage_mM)
xa_observed


# ── Distribution comparison ────────────────────────────────────────────────────

dists <- c("exponential", "weibull", "gompertz", "lognormal", "gengamma")

nested_models <- data_long2 |>
  nest() |>
  mutate(data = map(data, ~ mutate(.x, dosage_mM = fct_drop(dosage_mM)))) |>
  expand_grid(distribution = dists) |>
  mutate(
    model = map2(
      data,
      distribution,
      ~ fit_flexsurv(
        Surv(`Hours`, event) ~ experiment +
          dosage_mM +
          genotype +
          experiment:genotype,
        data = .x,
        dist = .y
      )
    ),
    aic = map_dbl(model, ~ if (!is.null(.x)) AIC(.x) else NA_real_)
  ) |>
  arrange(aic)

# ── Pooled model (as originally reported, now with the WT C- data) ────────────
# Arm (experiment) is a fixed blocking effect, so dosage effects are estimated
# within arms and assumed common to both arms and genotypes.

model_4 <- flexsurvreg(
  Surv(`Hours`, event) ~ experiment +
    dosage_mM +
    genotype,
  data = data_long2,
  dist = "weibull"
)

# ── Per-arm analysis of heterozygotes (primary; reviewer comment 4) ──────────
# Each treatment within an arm was a single cage, so cage (replicate) effects
# cannot be separated from treatment effects and no replicate term can be
# fitted. The two arms are independent experiments, so each is analysed
# separately and the XA effect is judged by its consistency across them.
# Heterozygotes only: across both arms just 1 of 112 WT females died (none in
# the sugar arm), so WT contribute no information on rescue and would cause
# complete separation; WT are shown by Kaplan-Meier curves only in Figure 4.

data_het <- data_long2 |>
  filter(genotype == "Het")

xa_arm_models <- data_het |>
  split(~experiment) |>
  map(\(d) {
    flexsurvreg(Surv(Hours, event) ~ dosage_mM, data = d, dist = "weibull")
  })

# Planned contrasts within each arm: XA vs NaOH vehicle (primary), XA vs
# untreated, and vehicle vs untreated (Holm-adjusted across the three)
xa_contrast_pairs <- tibble(
  level = c("6mM", "6mM", "C-"),
  ref = c("C-", "0mM", "0mM"),
  contrast = c("XA vs vehicle (C-)", "XA vs untreated", "C- vs untreated")
)

xa_arm_contrasts <- xa_arm_models |>
  imap(\(m, arm) {
    xa_contrast_pairs |>
      mutate(family = arm) |>
      contrast_levels(model = m, factor_name = "dosage_mM") |>
      mutate(experiment = arm, .before = 1)
  }) |>
  list_rbind()
xa_arm_contrasts

# Overall effect of treatment within each arm (likelihood-ratio test, 2 df)
xa_arm_lrt <- xa_arm_models |>
  imap(\(m, arm) {
    null <- flexsurvreg(
      Surv(Hours, event) ~ 1,
      data = filter(data_het, experiment == arm),
      dist = "weibull"
    )
    chisq <- 2 * (m$loglik - null$loglik)
    df <- m$npars - null$npars
    tibble(
      experiment = arm,
      chisq = chisq,
      df = df,
      p = pchisq(chisq, df, lower.tail = FALSE)
    )
  }) |>
  list_rbind()
xa_arm_lrt

# Survival probability at 72 h by arm and treatment (heterozygotes)
xa_surv_72 <- xa_arm_models |>
  imap(\(m, arm) {
    summary(
      m,
      newdata = tibble(dosage_mM = levels(data_het$dosage_mM)),
      type = "survival",
      t = 72,
      tidy = TRUE
    ) |>
      as_tibble() |>
      mutate(experiment = arm, .before = 1)
  }) |>
  list_rbind()
xa_surv_72

# ── Combined evidence across arms (supporting) ────────────────────────────────
# Heterozygote-only pooled model with arm as a blocking factor, and a test of
# whether the treatment effect differs between arms (arm x treatment).

xa_pooled_het <- flexsurvreg(
  Surv(Hours, event) ~ experiment + dosage_mM,
  data = data_het,
  dist = "weibull"
)
xa_pooled_het_int <- flexsurvreg(
  Surv(Hours, event) ~ experiment * dosage_mM,
  data = data_het,
  dist = "weibull"
)

xa_interaction_lrt <- tibble(
  chisq = 2 * (xa_pooled_het_int$loglik - xa_pooled_het$loglik),
  df = xa_pooled_het_int$npars - xa_pooled_het$npars
) |>
  mutate(p = pchisq(chisq, df, lower.tail = FALSE))
xa_interaction_lrt

xa_pooled_contrasts <- xa_contrast_pairs |>
  mutate(family = "pooled") |>
  contrast_levels(model = xa_pooled_het, factor_name = "dosage_mM")
xa_pooled_contrasts


# ── Predicted survival curves (per-arm heterozygote models) ──────────────────

surv_df <- xa_arm_models |>
  imap(\(m, arm) {
    pred_surv(
      m,
      newdata = tibble(dosage_mM = levels(data_het$dosage_mM)),
      t_max = max(filter(data_het, experiment == arm)$Hours)
    ) |>
      mutate(experiment = arm, genotype = "Het")
  }) |>
  list_rbind() |>
  mutate(
    dosage_mM = factor(dosage_mM, levels = levels(data_long2$dosage_mM)),
    genotype = factor(genotype, levels = c("WT", "Het"))
  )

km_model_4_df <- tidy_km(
  data_long2,
  Surv(Hours, event) ~ experiment + dosage_mM + genotype
) |>
  separate_wider_delim(
    strata,
    delim = ",",
    names = c("experiment", "dosage_mM", "genotype"),
    too_many = "merge"
  ) |>
  mutate(across(
    .cols = c("experiment", "dosage_mM", "genotype"),
    .fns = str_trim
  )) |>
  mutate(
    dosage_mM = factor(dosage_mM, levels = levels(data_long2$dosage_mM)),
    genotype = factor(genotype, levels = c("WT", "Het"))
  )


# ── Shared plot scales and theme ──────────────────────────────────────────────

dose_colours <- c(
  "0mM" = "lightgrey",
  "C-" = "#56B4E9",
  "6mM" = "darkgreen"
)

suppl_scales <- list(
  scale_colour_manual(
    values = dose_colours,
    breaks = names(dose_colours),
    name = "Treatment"
  ),
  scale_fill_manual(values = dose_colours, name = "Treatment"),
  estimate_linetype_layers(),
  scale_y_continuous(limits = c(0, 1)),
  labs(x = "Time (Hours)", y = "Survival Probability"),
  theme_bw(base_size = 12),
  theme(
    axis.text.x = element_text(),
    panel.grid.minor = element_blank(),
    legend.key = element_blank()
  )
)

# ── Figure 4: facet by genotype x supplementation route ──────────────────────

XA_suppl <- ggplot(
  surv_df,
  aes(
    x = time,
    y = est,
    colour = dosage_mM,
    fill = dosage_mM
  )
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, colour = NA) +
  geom_line(aes(linetype = lt_model), linewidth = 0.8) +
  geom_step(
    data = km_model_4_df,
    aes(x = time, y = estimate, colour = dosage_mM, linetype = lt_km),
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  facet_grid(genotype ~ experiment) +
  suppl_scales

XA_suppl
