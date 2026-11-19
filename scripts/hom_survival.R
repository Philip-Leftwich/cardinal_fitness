source(here::here("scripts", "packages.R"))

path <- here::here("data", "SurvivalData.xlsx")

eclosion_df <- path |>
  excel_sheets() |>
  set_names() |>
  map(\(s) read_excel(path = path, sheet = s)) |>
  list_rbind(names_to = "line")

line_order <- c("QA383P", "1759")

eclosion_df <- eclosion_df |>
  mutate(
    Genotype = case_when(
      Genotype == "SDA-500" ~ "WT",
      Genotype == "CdKO Het" ~ "HET",
      Genotype == "CdKO Hom" ~ "HOM",
      .default = as.character(Genotype)
    )
  ) |>
  mutate(Genotype = factor(Genotype, levels = c("HET", "WT", "HOM"))) |>
  mutate(Sex = str_to_lower(Sex)) |>
  mutate(line = factor(line, levels = line_order))

n_total <- 200  # total pupae scored per vial (binomial denominator)

# Biological plausibility checks ----
# Hard checks: adult counts must be non-negative and cannot exceed the total
# number of pupae scored; verified before modelling (Supplementary figure)
eclosion_df <- eclosion_df |>
  verify(Adults >= 0) |>                   # eclosion count cannot be negative
  verify(Adults <= n_total)                 # eclosion cannot exceed total pupae scored


model <- glm(
  cbind(`Adults`, n_total - `Adults`) ~ line * Genotype * Sex,
  family = binomial,
  data = eclosion_df
)
summary(model)

# Observed counts. Each tray held 200 mixed-sex larvae and adults were scored
# by sex, so each sex-specific proportion has a maximum of ~50% under a 1:1
# sex ratio. Note the WT group differs by line: cd[g384] = non-transgenic
# siblings; cd[384R] = SDA-500.
eclosion_observed <- eclosion_df |>
  summarise(
    trays = n(),
    larvae = sum(Larvae),
    pupae = sum(`Pupae...5`),
    adults = sum(Adults),
    .by = c(line, Genotype, Sex)
  ) |>
  mutate(prop_adults = adults / larvae)
eclosion_observed

# Each genotype versus WT within line and sex (odds ratios; Dunnett-adjusted
# p-values, unadjusted 95% CIs)
eclosion_emm <- emmeans::emmeans(model, ~ Genotype | line + Sex)
eclosion_vs_wt <- emmeans::contrast(
  eclosion_emm,
  method = "trt.vs.ctrl",
  ref = which(levels(eclosion_df$Genotype) == "WT")
) |>
  summary(infer = c(FALSE, TRUE)) |>
  as_tibble() |>
  mutate(
    ratio = exp(estimate),
    ratio.LCL = exp(estimate - qnorm(0.975) * SE),
    ratio.UCL = exp(estimate + qnorm(0.975) * SE)
  )
eclosion_vs_wt
