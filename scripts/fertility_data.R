source(here::here("scripts", "packages.R"))
source(here::here("scripts", "functions", "fertility_data_functions.R"))

# Read survival data ====

filepath <- here::here("data", "FertilityData.xlsx")

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

# Biological plausibility checks ----
# Hard checks: verify absolute range constraints and cross-column logic
# before any downstream modelling (produces Figure 1)
fertility_data <- fertility_data |>
  verify(eggs >= 0 | is.na(eggs)) |>           # egg counts cannot be negative
  verify(larvae >= 0 | is.na(larvae)) |>        # larval counts cannot be negative
  verify(is.na(eggs) | is.na(larvae) | larvae <= eggs)  # larvae cannot exceed eggs

# Cross label order
cross_order <- c(
  "WT",
  "Female \nHET",
  "Male \nHET",
  "Female \nHOM",
  "Male \nHOM"
)

line_order <- c("QA383P", "2360B5", "1759")

fertility_data <- fertility_data |>
  mutate(
    cross = case_when(
      # 2360B5
      line == "2360B5" & cross == "HETxSDA" ~ "Female \nHET",
      line == "2360B5" & cross == "SDAxHET" ~ "Male \nHET",
      line == "2360B5" & cross == "HOMxSDA" ~ "Female \nHOM",
      line == "2360B5" & cross == "SDAxHOM" ~ "Male \nHOM",
      # QA383P  cross labels
      line == "QA383P" & cross == "CdKO(het)xSDA" ~ "Male \nHET",
      line == "QA383P" & cross == "SDAxCdKO(het)" ~ "Female \nHET",
      line == "QA383P" & cross == "CdKO(hom)xSDA" ~ "Male \nHOM",
      line == "QA383P" & cross == "SDAxCdKO(hom)" ~ "Female \nHOM",
      # 1759 cross labels
      line == "1759" & cross == "HETxSDA" ~ "Male \nHET",
      line == "1759" & cross == "SDAxHET" ~ "Female \nHET",
      line == "1759" & cross == "HOMxSDA" ~ "Male \nHOM",
      line == "1759" & cross == "SDAxHOM" ~ "Female \nHOM",
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
  y_label = "Egg count"
)

p_fecundity

# Fertility ====

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
      c("1759", "QA383P") |
      line == "2360B5" & cross == "Male \nHET" |
      line == "2360B5" & cross == "WT"
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
  y_scale = scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  )
)

p_fertility
