source("scripts/packages.R")

path <- "clean_data/SurvivalData.xlsx"

data <- path |>
  excel_sheets() |>
  set_names() |>
  map_df(~ read_excel(path = path, sheet = .x), .id = "line")

line_order <- c("QA383P", "1759")

data <- data |>
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


model <- glm(
  cbind(`Adults`, 200 - `Adults`) ~ line * Genotype * Sex,
  family = binomial,
  data = data
)
summary(model)
