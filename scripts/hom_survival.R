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
