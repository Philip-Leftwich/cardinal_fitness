source(here::here("scripts", "packages.R"))
source(here::here("scripts", "smurf_functions.R"))

# Load data ====

path <- "data/smurf.xlsx"

line_labels <- c(
  "KI" = "italic(cd)^g225",
  "KO" = "italic(cd)^{225*R}",
  "WT" = "wildtype"
)


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
  ) |>
  mutate(asymp.UCL = if_else(asymp.UCL > 0.990, 0, asymp.UCL)) # cap CI

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
  scale_colour_manual(
    values = c(
      "lightgray",
      "#8BABD3",
      "#FFA040"
    ),
    name = "Genotype"
  ) +
  labs(
    x = "Genotype",
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
  scale_colour_manual(
    values = c(
      "lightgray",
      "#8BABD3",
      "#FFA040"
    ),
    name = "Genotype"
  ) +
  scale_fill_manual(
    values = c(
      "lightgray",
      "#8BABD3",
      "#FFA040"
    ),
    name = "Genotype"
  ) +
  labs(
    x = "Genotype",
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
  filter(delivery == "sucrose") |>
  mutate(line = str_to_upper(line))

# Weibull model + predicted survival curves
weibull_after <- fit_weibull_surv(
  data = survival_after_smurf,
  formula = Surv(hours, survival) ~ line,
  t_max = 120
)
model_smurf_surv <- weibull_after$model
surv_df_smurf <- weibull_after$surv_df

surv_df_smurf <- surv_df_smurf |>
  mutate(lcl = if_else(est > 0.999, 1, lcl)) |>
  mutate(line = factor(line, levels = genotype_order))


# Kaplan-Meier estimates for overlay
km_df_smurf <- tidy_km(survival_after_smurf, Surv(hours, survival) ~ line) |>
  rename(line = strata) |>
  mutate(line = factor(line, levels = genotype_order))

# Visualisation
survival_after_smurf_plot <- plot_surv(
  surv_df = surv_df_smurf,
  km_df = km_df_smurf,
  colour_var = "line",
  colour_scale = scale_colour_manual(
    values = c(
      "lightgray",
      "#8BABD3",
      "#FFA040"
    ),
    name = "Genotype"
  ),
  fill_scale = scale_fill_manual(
    values = c(
      "lightgray",
      "#8BABD3",
      "#FFA040"
    ),
    name = "Genotype"
  ),
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
smurf_colours <- c("Smurf-positive" = "#4B779D", "Smurf-negative" = "#606060")

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
