source(here::here("scripts", "packages.R"))
source(here::here("scripts", "fertility_data_functions.R"))

# Read survival data ====

filepath <- here::here("clean_data", "FertilityData.xlsx")

sheets <- excel_sheets(filepath)

selected_sheets <- sheets[c(1, 2, 4)]

data_list <- map(
  selected_sheets,
  ~ read_excel(
    filepath,
    sheet = .x,
    na = c("ND", "N/A", "_", "-"),
    .name_repair = "unique"
  ) |>
    janitor::clean_names() |>
    mutate(line = paste(.x))
)


# Standardise sheets
data_list[[1]] <- data_list[[1]] |>
  unite("plate_well", c(plate, well), sep = "+") |>
  select(plate_well, cross, rep, eggs, larvae, line)

data_list[[2]] <- data_list[[2]] |>
  rename("rep" = "replicate", "eggs" = "egg_num", "larvae" = "larvae_num") |>
  select(plate_well, cross, rep, eggs, larvae, line)

data_list[[3]] <- data_list[[3]] |>
  mutate(plate = cumsum(well_position == "A1")) |>
  unite("plate_well", c(plate, well_position), sep = "+") |>
  rename(
    "rep" = "replicate_cage",
    "eggs" = "egg_num",
    "larvae" = "larvae_num"
  ) |>
  select(plate_well, cross, rep, eggs, larvae, line)

# recombine data
fertility_data <- rbind(data_list[[1]], data_list[[2]], data_list[[3]])

# Cross label order
cross_order <- c("WT", "HETxSDA", "SDAxHET", "HOMxSDA", "SDAxHOM")

line_order <- c("QA383P", "2360B5", "1759")

fertility_data <- fertility_data |>
  mutate(
    cross = case_when(
      cross == "CdKO(het)xSDA" ~ "HETxSDA",
      cross == "SDAxCdKO(het)" ~ "SDAxHET",
      cross == "CdKO(hom)xSDA" ~ "HOMxSDA",
      cross == "SDAxCdKO(hom)" ~ "SDAxHOM",
      cross %in% c("WTxSDA", "SDAxWT", "SDAxSDA") ~ "WT",
      .default = as.character(cross)
    )
  ) |>
  drop_na(cross) |>
  mutate(cross = factor(cross, levels = cross_order)) |>
  mutate(line = factor(line, levels = line_order)) |>
  unite("line_cross", line, cross, remove = FALSE)

egg_model <- glmmTMB(
  eggs ~ line * cross + (1 | rep / plate_well),
  family = nbinom1(),
  ziformula = ~1,
  data = fertility_data
)

# Get emmeans predictions on response scale
egg_emm <- get_emm(egg_model, cross_order)

# Raw egg counts (exclude NA eggs)
raw_eggs <- fertility_data |>
  filter(!is.na(eggs)) |>
  mutate(cross = factor(cross, levels = cross_order))

# Plot fecundity predictions
p_fecundity <- plot_fertility(
  emm_data = egg_emm,
  raw_data = raw_eggs,
  raw_y = eggs,
  y_label = "Egg count",
  plot_title = "Fecundity: number of eggs laid"
)

p_fecundity

#Fertility ====

# Almost all 2360 crosses produced no eggs so hatching model fails to converge
fertility_data_no2360 <- fertility_data |> filter(line != "2360B5")

larvae_model <- glmmTMB(
  cbind(larvae, (eggs - larvae)) ~ line * cross + (1 | line / rep / plate_well),
  family = binomial,
  data = fertility_data_no2360
)


# Create a combined cross and line category to analyse the 2360 group
fertility_data_combined <- fertility_data |>
  filter(
    line %in%
      c("1759", "QA383P", "WT") |
      (line == "2360B5" & cross == "SDAxHET")
  ) |>
  unite("line_cross", line, cross, remove = FALSE)

larvae_model2 <- glmmTMB(
  cbind(larvae, (eggs - larvae)) ~ line * cross + (1 | line / rep / plate_well),
  family = binomial,
  data = fertility_data_combined
)

# Plot fertility predictions ====

# Get emmeans predictions on response scale
larvae_emm <- get_emm(larvae_model2, cross_order)


# Raw proportions (exclude zero-egg wells)
raw_props <- fertility_data_combined |>
  filter(!is.na(larvae), eggs > 0) |>
  mutate(
    prop = larvae / eggs,
    cross = factor(cross, levels = cross_order)
  )

# Plot
p_fertility <- plot_fertility(
  emm_data = larvae_emm,
  raw_data = raw_props,
  raw_y = prop,
  y_label = "Hatching rate (larvae / eggs)",
  plot_title = "Fertility: proportion of eggs hatching to larvae",
  y_scale = scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  )
)

p_fertility
