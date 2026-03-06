library(tidyverse)
library(here)
library(janitor)
library(readxl)
library(patchwork)
library(emmeans)
library(performance)
library(glmmTMB)
library(sjPlot)

# Read survival data ====

filepath <- here::here("clean_data", "FertilityData.xlsx")

sheets <- excel_sheets(filepath)

selected_sheets <- sheets[c(1,2,4)]


## Data sheets are inconsistently formatted 
# read in as a list
data_list <- map(selected_sheets, 
                     ~ read_excel(filepath, 
                                  sheet = .x,
                                  na = c("ND", "N/A", "_", "-"),
                                  .name_repair = make_clean_names) %>% 
                       mutate(line = paste(.x)))

# Standardise sheets
data_list[[1]] <- data_list[[1]] %>% 
  unite("plate_well", c(plate, well), sep = "+") %>% 
  select(plate_well, cross, rep, eggs, larvae, line)

data_list[[2]] <- data_list[[2]] %>% 
  rename("rep" = "replicate",
         "eggs" = "egg_num",
         "larvae" = "larvae_num") %>% 
  select(plate_well, cross, rep, eggs, larvae, line)

data_list[[3]] <- data_list[[3]] %>% 
  mutate(plate = cumsum(well_position == "A1")) %>% 
  unite("plate_well", c(plate, well_position), sep = "+") %>% 
  rename("rep" = "replicate_cage",
         "eggs" = "egg_num",
         "larvae" = "larvae_num") %>% 
  select(plate_well, cross, rep, eggs, larvae, line)

# recombine data
fertility_data <- rbind(data_list[[1]], data_list[[2]], data_list[[3]])

# Cross label order
cross_order <- c("WT","HETxSDA","SDAxHET", "HOMxSDA", "SDAxHOM")

line_order <- c("QA383PB5", "2360B5", "1759B5")

fertility_data <- fertility_data %>% 
  mutate(cross = case_when(
    cross == "CdKO(het)xSDA" ~ "HETxSDA",
    cross == "SDAxCdKO(het)" ~ "SDAxHET", 
    cross == "CdKO(hom)xSDA" ~ "HOMxSDA", 
    cross == "SDAxCdKO(hom)" ~ "SDAxHOM", 
    cross %in% c("WTxSDA", "SDAxWT", "SDAxSDA") ~ "WT",
    .default = as.character(cross))) %>% 
    drop_na(cross) |> 
  mutate(cross = factor(cross, levels = cross_order)) |> 
  mutate(line = factor(line, levels = line_order)) |>
  unite("line_cross", line,cross,
remove = FALSE)
# remove wildtype crosses
#fertility_data <- fertility_data %>% 
#  filter(!cross %in% c("WTxSDA", "SDAxWT", "SDAxSDA"))

# Fecundity ====



egg_model <- glmmTMB(eggs ~ line * cross + (1|line/rep/plate_well), 
                     family = nbinom1(), 
                     ziformula = ~ 1,
                     data = fertility_data)

# Get emmeans predictions on response scale
egg_emm <- emmeans::emmeans(egg_model, ~ line + cross,
                            type = "response") 
  as_tibble()



egg_emm <- egg_emm |>
  mutate(cross = factor(cross, levels = cross_order)) |> 
  drop_na()

# Raw egg counts (exclude NA eggs)
raw_eggs <- fertility_data |>
  filter(!is.na(eggs)) |>
  mutate(cross = factor(cross, levels = cross_order))

# Plot fecundity predictions
p_fecundity <- ggplot(egg_emm, aes(x = cross, y = response, colour = line, group = line)) +
  # Raw data jittered points
  geom_point(
    data = raw_eggs,
    aes(y = eggs, colour = line),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    alpha = 0.3, size = 1.5, shape = 16
  ) +
  # Model predictions
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  # 95% CI
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.6),
    width = 0.25, linewidth = 0.8
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Line") +
  labs(
    x = "Cross",
    y = "Egg count",
    title = "Fecundity: number of eggs laid"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_fecundity

#Fertility ====

# Almost all 2360 crosses produced no eggs so hatching model fails to converge
fertility_data_no2360 <-  fertility_data %>% filter(line != "2360B5")

larvae_model <- glmmTMB(cbind(larvae,(eggs-larvae)) ~ line * cross + (1|line/rep/plate_well),
                        family = binomial,
                        data = fertility_data_no2360)


# Create a combined cross and line category to analyse the 2360 group
fertility_data_combined <- fertility_data %>% 
  filter(line %in% c("1759", "QA383P", "WT") | (line == "2360B5" & cross == "SDAxHET")) %>% 
  unite("line_cross", line,cross,
remove = FALSE)

larvae_model2 <- glmmTMB(cbind(larvae,(eggs-larvae)) ~ line * cross + (1|line/rep/plate_well),
                        family = binomial,
                        data = fertility_data_combined)

# Plot fertility predictions ====

# Get emmeans predictions on response scale
larvae_emm <- emmeans::emmeans(larvae_model2, ~ line + cross, type = "response") |>
  as.data.frame()
 # separate(line_cross, into = c("line", "cross"), sep = "_")

# Raw proportions (exclude zero-egg wells)
raw_props <- fertility_data_combined |>
  filter(!is.na(larvae), eggs > 0) |>
  mutate(prop = larvae / eggs)

# Cross label order: HOM crosses last to match biological expectation
#cross_order <- c("HETxSDA", "SDAxHET", "HOMxSDA", "SDAxHOM")

larvae_emm <- larvae_emm |>
  mutate(cross = factor(cross, levels = cross_order))

raw_props <- raw_props |>
  mutate(cross = factor(cross, levels = cross_order))

# Plot
p_fertility <- ggplot(larvae_emm, aes(x = cross, y = response, colour = line, group = line)) +
  # Raw data jittered points
  geom_point(
    data = raw_props,
    aes(y = prop, colour = line),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    alpha = 0.3, size = 1.5, shape = 16
  ) +
  # Model predictions
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  # 95% CI
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.6),
    width = 0.25, linewidth = 0.8
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_brewer(palette = "Dark2", name = "Line") +
  labs(
    x = "Cross",
    y = "Hatching rate (larvae / eggs)",
    title = "Fertility: proportion of eggs hatching to larvae"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_fertility
