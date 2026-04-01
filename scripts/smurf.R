source(here::here("scripts", "packages.R"))

# Load data ====

path <- "data/smurf.xlsx"

line_labels <- c(
  "ki" = "italic(cd)^g225",
  "ko" = "italic(cd)^225",
  "wt" = "wildtype"
)

# ── Shared helper functions ────────────────────────────────────────────────────

# Fit Weibull model and return tidy predicted survival curve data frame
fit_weibull_surv <- function(data, formula, t_max, n_t = 108) {
  model <- flexsurvreg(formula, data = data, dist = "weibull")
  pred_vars <- all.vars(formula)[-c(1, 2)] # predictor names from formula
  pred_data <- map(pred_vars, ~ unique(data[[.x]])) |>
    setNames(pred_vars) |>
    do.call(what = expand_grid, args = _)
  surv_pred <- summary(
    model,
    newdata = pred_data,
    type = "survival",
    t = seq(0, t_max, length.out = n_t),
    tidy = TRUE
  )
  list(model = model, surv_df = as_tibble(surv_pred))
}

# Tidy Kaplan-Meier fit, stripping strata prefix, with time-zero row added
tidy_km <- function(data, formula) {
  km_fit <- survfit(formula, data = data)
  strata_prefix <- paste0(all.vars(formula[[3]]), "=")
  tidy_df <- broom::tidy(km_fit) |>
    mutate(across(
      strata,
      ~ str_remove_all(.x, paste(strata_prefix, collapse = "|"))
    ))
  # Prepend time = 0, estimate = 1 for each stratum so step lines start at 1
  time_zero <- tidy_df |>
    distinct(strata) |>
    mutate(time = 0, estimate = 1)
  bind_rows(time_zero, tidy_df) |>
    arrange(strata, time)
}

# Build Weibull + KM survival plot
plot_surv <- function(
  surv_df,
  km_df,
  colour_var,
  colour_scale,
  fill_scale,
  extra_layers = NULL
) {
  p <- ggplot(surv_df, aes(x = time, y = est, fill = .data[[colour_var]])) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
    geom_line(aes(colour = .data[[colour_var]])) +
    geom_step(
      data = km_df,
      aes(x = time, y = estimate, colour = .data[[colour_var]]),
      linetype = "dashed",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 120, by = 24)) +
    colour_scale +
    fill_scale +
    labs(
      x = "Time (Hours)",
      y = "Survival Probability",
      caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank()
    )
  if (!is.null(extra_layers)) {
    p <- p + extra_layers
  }
  p
}


data <- read_excel(
  path = path,
  sheet = "Smurf screening",
  .name_repair = janitor::make_clean_names,
  skip = 1
) |>
  rename("id" = "x") |>
  select(!contains("x")) |>
  rename_with(
    .cols = 2:25,
    .fn = ~ case_when(
      str_detect(., "_2") ~ str_replace(., "_2", "_KO_TRP"),
      str_detect(., "_3") ~ str_replace(., "_3", "_KO_XA"),
      str_detect(., "_4") ~ str_replace(., "_4", "_KO_NAOH"),
      str_detect(., "_5") ~ str_replace(., "_5", "_KI_SUCROSE"),
      str_detect(., "_6") ~ str_replace(., "_6", "_KI_TRP"),
      str_detect(., "_7") ~ str_replace(., "_7", "_KI_XA"),
      str_detect(., "_8") ~ str_replace(., "_8", "_KI_NAOH"),
      str_detect(., "_9") ~ str_replace(., "_9", "_WT_SUCROSE"),
      str_detect(., "_10") ~ str_replace(., "_10", "_WT_TRP"),
      str_detect(., "_11") ~ str_replace(., "_11", "_WT_XA"),
      str_detect(., "_12") ~ str_replace(., "_12", "_WT_NAOH"),
      TRUE ~ str_c(., "_KO_SUCROSE")
    )
  ) |>
  rename_with(
    .cols = 2:25,
    .fn = ~ str_replace(., "smurf_", "")
  ) |>
  pivot_longer(cols = 2:25, names_to = "condition", values_to = "values") |>
  separate(condition, c("smurf", "genotype", "delivery")) |>
  mutate(id = str_remove_all(id, "[0-9]")) |>
  drop_na(id)


data_wide <- data |>
  group_by(id, smurf, genotype, delivery) |>
  summarise(n = sum(values, na.rm = T)) |>
  pivot_wider(names_from = smurf, values_from = n) |>
  ungroup() |>
  mutate(delivery = fct_relevel(delivery, "SUCROSE")) |>
  mutate(genotype = fct_relevel(genotype, "WT"))

data_wide <- data_wide |>
  filter(delivery %in% "SUCROSE")

# Model=====

smurf_model <- glm(
  cbind(positive, negative) ~ genotype,
  data = data_wide,
  family = binomial(link = "logit")
)

drop1(smurf_model, test = "Chisq")

summary(smurf_model)

# Plot smurf predictions ====

# Get emmeans predictions on response scale
smurf_emm <- emmeans::emmeans(smurf_model, ~genotype, type = "response") |>
  as.data.frame()

# Factor levels

genotype_order <- c("WT", "KI", "KO")

smurf_emm <- smurf_emm |>
  mutate(
    genotype = factor(genotype, levels = genotype_order)
  )

# Raw proportions
raw_smurf <- data_wide |>
  mutate(
    prop = positive / (positive + negative),
    genotype = factor(genotype, levels = genotype_order)
  )

# Plot
p_smurf <- ggplot(
  smurf_emm,
  aes(x = genotype, y = prob, colour = genotype, group = genotype)
) +

  # 95% CI
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_nudge(x = .1),
    width = 0.15,
    linewidth = 0.8
  ) +
  # Model predictions
  geom_point(
    position = position_nudge(x = .1),
    size = 3
  ) +
  # Raw data
  geom_point(
    data = raw_smurf,
    aes(x = genotype, y = prop, colour = genotype),
    position = position_jitter(width = .15),
    alpha = 0.3,
    size = 1.5,
    shape = 16
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_brewer(palette = "Dark2", name = "Genotype") +
  labs(
    x = "Delivery",
    y = "Smurf rate",
    title = "Smurf assay: proportion of smurf-positive flies"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_smurf


p_smurf2 <- ggplot(
  smurf_emm,
  aes(
    x = genotype,
    y = prob,
    colour = genotype,
    fill = genotype,
    group = genotype
  )
) +
  # Model predictions
  geom_col(alpha = .4) +

  # 95% CI
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.15,
    linewidth = 0.8
  ) +
  # Raw data
  geom_point(
    data = raw_smurf,
    aes(x = genotype, y = prop, colour = genotype),
    position = position_jitter(width = .15),
    alpha = 0.3,
    size = 1.5,
    shape = 16
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_brewer(palette = "Dark2", name = "Genotype") +
  scale_fill_brewer(palette = "Dark2", name = "Genotype") +
  labs(
    x = "Delivery",
    y = "Smurf rate",
    title = "Smurf assay: proportion of smurf-positive flies"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_smurf2


#########################
# Survival after Smurf assay ====
#########################

path <- "clean_data/SmurfAssayData.xlsx"
sheets <- c("Survival")

survival_after_smurf <- read_excel(
  path = path,
  sheet = sheets,
  .name_repair = janitor::make_clean_names
)

survival_after_smurf <- survival_after_smurf |>
  pivot_longer(
    cols = -c(x, hours),
    names_to = "treatment",
    values_to = "survival"
  ) |>
  drop_na(survival) |>
  separate(treatment, into = c("line", "delivery"), sep = "_") |>
  filter(delivery == "sucrose")

# Weibull model + predicted survival curves
weibull_after <- fit_weibull_surv(
  data = survival_after_smurf,
  formula = Surv(hours, survival) ~ line,
  t_max = 120
)
model_smurf_surv <- weibull_after$model
surv_df_smurf <- weibull_after$surv_df

# Kaplan-Meier estimates for overlay
km_df_smurf <- tidy_km(survival_after_smurf, Surv(hours, survival) ~ line) |>
  rename(line = strata)

# Visualisation
survival_after_smurf_plot <- plot_surv(
  surv_df = surv_df_smurf,
  km_df = km_df_smurf,
  colour_var = "line",
  colour_scale = scale_colour_brewer(palette = "Dark2", name = "Genotype"),
  fill_scale = scale_fill_brewer(palette = "Dark2", name = "Genotype"),
  extra_layers = facet_wrap(
    ~line,
    labeller = as_labeller(line_labels, label_parsed)
  )
)

survival_after_smurf_plot


#########################
# Survival by smurf type ====
#########################

sheets <- c("SurvivalSmurf")

survival_by_smurf <- read_excel(
  path = path,
  sheet = sheets,
  .name_repair = janitor::make_clean_names
)

survival_by_smurf <- survival_by_smurf |>
  pivot_longer(
    cols = -c(x, hours),
    names_to = "treatment",
    values_to = "survival"
  ) |>
  mutate(
    smurf = case_when(
      treatment == "ki_sucrose" ~ "Smurf-positive",
      treatment == "ki_sucrose_2" ~ "Smurf-negative",
      TRUE ~ NA_character_
    )
  ) |>
  drop_na(survival, smurf) # drop both NA survival AND NA smurf


# Weibull model + predicted survival curves
smurf_colours <- c("Smurf-positive" = "#E69F00", "Smurf-negative" = "#4D4D4D")

weibull_by <- fit_weibull_surv(
  data = survival_by_smurf,
  formula = Surv(hours, survival) ~ smurf,
  t_max = 120
)
model_smurf_by <- weibull_by$model
surv_df_by_smurf <- weibull_by$surv_df

# Kaplan-Meier estimates for overlay
km_df_by_smurf <- tidy_km(survival_by_smurf, Surv(hours, survival) ~ smurf) |>
  rename(smurf = strata)

# Visualisation
survival_by_smurf_plot <- plot_surv(
  surv_df = surv_df_by_smurf,
  km_df = km_df_by_smurf,
  colour_var = "smurf",
  colour_scale = scale_colour_manual(
    values = smurf_colours,
    name = "Smurf status"
  ),
  fill_scale = scale_fill_manual(values = smurf_colours, name = "Smurf status")
)

survival_by_smurf_plot
