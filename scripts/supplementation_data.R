source(here::here("scripts", "packages.R"))
source(here::here("scripts", "supplementation_data_functions.R"))

# Load data ====

path <- "clean_data/XASupplementationData.xlsx"

data <- path |>
  excel_sheets() |>
  set_names() |>
  map_df(~ read_excel(path = path, sheet = .x), .id = "experiment") |>

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

data_long2 <- data_long |>
  separate(treatment, into = c("genotype", "dosage"), sep = "\\(") |>
  mutate(genotype = str_trim(genotype)) |>
  mutate(dosage = str_remove(dosage, "\\)")) |>
  mutate(experiment = str_remove(experiment, "XA ")) |>
  separate(dosage, c("dosage", "c-")) |>
  mutate("c_minus" = factor(if_else(`c-` == "C", "-", "+", missing = "+"))) |>
  mutate(
    dosage = factor(str_remove(dosage, "mM"), levels = c("0", "6", "24"))
  ) |>
  rename("dosage_mM" = dosage)


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

# ── Weibull model selection ────────────────────────────────────────────────────
# Models are progressively simplified — model4 selected as best fit

weibull_formulas <- list(
  model = Surv(`Hours`, event) ~ experiment +
    dosage_mM +
    genotype +
    c_minus +
    experiment:dosage_mM +
    experiment:genotype +
    experiment:c_minus +
    dosage_mM:genotype +
    genotype:c_minus,
  model2 = Surv(`Hours`, event) ~ experiment +
    dosage_mM +
    genotype +
    c_minus +
    experiment:dosage_mM +
    experiment:genotype +
    experiment:c_minus,
  model3 = Surv(`Hours`, event) ~ experiment +
    dosage_mM +
    genotype +
    c_minus +
    experiment:dosage_mM +
    experiment:c_minus,
  model4 = Surv(`Hours`, event) ~ experiment +
    dosage_mM +
    genotype +
    c_minus +
    experiment:c_minus
)

model_4 <- flexsurvreg(
  Surv(`Hours`, event) ~ experiment +
    dosage_mM +
    genotype +
    c_minus +
    experiment:c_minus,
  data = data_long2,
  dist = "weibull"
)


# ── Predicted survival curves ─────────────────────────────────────────────────

pred_data <- expand_grid(
  experiment = unique(data_long2$experiment),
  dosage_mM = factor(c(0, 6, 24), levels = c("0", "6", "24")),
  genotype = unique(data_long2$genotype),
  c_minus = unique(data_long2$c_minus)
)

surv_df <- pred_surv(
  model_4,
  newdata = pred_data,
  t_max = max(data_long2$Hours)
)

# ── Shared plot scales and theme ──────────────────────────────────────────────

dose_colours <- c("0" = "#00C9A7", "6" = "#9B72CF", "24" = "#5C2D91")

#dose_colours <- setNames(RColorBrewer::brewer.pal(3, "Purples"), c("0", "6", "24"))

suppl_scales <- list(
  scale_colour_manual(values = dose_colours, name = "Dosage (mM)"),
  scale_fill_manual(values = dose_colours, name = "Dosage (mM)"),
  scale_linetype_manual(values = c("+" = "solid", "-" = "dashed"), name = "C-"),
  scale_y_continuous(limits = c(0, 1)),
  labs(x = "Time (Hours)", y = "Survival Probability"),
  guides(
    linetype = guide_legend(
      override.aes = list(fill = NA, colour = "black")
    )
  ),
  theme_bw(base_size = 12),
  theme(
    axis.text.x = element_text(),
    panel.grid.minor = element_blank(),
    legend.key = element_blank()
  )
)

# ── Plot 1: linetype for c_minus, facet by experiment × genotype ──────────────

XA_suppl_v1 <- ggplot(
  surv_df,
  aes(
    x = time,
    y = est,
    colour = dosage_mM,
    fill = dosage_mM,
    linetype = c_minus,
    group = interaction(dosage_mM, c_minus)
  )
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(experiment ~ genotype) +
  suppl_scales

XA_suppl_v1

# ── Plot 2: ghost lines — c_minus comparison within each panel ────────────────

surv_df_ghost <- bind_rows(
  surv_df |> mutate(panel_c = c_minus, alpha_val = 1),
  surv_df |>
    mutate(panel_c = if_else(c_minus == "+", "-", "+"), alpha_val = 0.25)
)

XA_suppl <- ggplot(
  surv_df_ghost,
  aes(
    x = time,
    y = est,
    colour = dosage_mM,
    fill = dosage_mM,
    linetype = c_minus,
    group = interaction(dosage_mM, c_minus)
  )
) +
  geom_ribbon(
    aes(ymin = lcl, ymax = ucl, alpha = I(alpha_val * 0.08)),
    colour = NA
  ) +
  geom_line(aes(alpha = I(alpha_val)), linewidth = 0.8) +
  facet_grid(experiment + panel_c ~ genotype) +
  suppl_scales +
  labs(
    caption = "Solid: focal C± group  |  Faded: opposite group for comparison"
  )

XA_suppl


genotype_colours <- c("WT" = "#4D4D4D", "Het" = "#E69F00")

XA_suppl_v3 <- ggplot(
  surv_df,
  aes(
    x = time,
    y = est,
    colour = genotype,
    fill = genotype,
    linetype = dosage_mM,
    group = interaction(genotype, dosage_mM)
  )
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    experiment ~ c_minus,
    labeller = labeller(c_minus = c("+" = "C+", "-" = "C-"))
  ) +
  scale_colour_manual(values = genotype_colours, name = "Genotype") +
  scale_fill_manual(values = genotype_colours, name = "Genotype") +
  scale_linetype_manual(
    values = c("0" = "dotted", "6" = "dashed", "24" = "solid"),
    name = "Dosage (mM)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time (Hours)", y = "Survival Probability") +
  guides(
    colour = guide_legend(override.aes = list(fill = NA)),
    linetype = guide_legend(override.aes = list(fill = NA, colour = "black"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.key = element_blank()
  )

XA_suppl_v3
