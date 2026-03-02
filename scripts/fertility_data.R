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

fertility_data <- fertility_data %>% 
  mutate(cross = case_when(
    cross == "CdKO(het)xSDA" ~ "HETxSDA",
    cross == "SDAxCdKO(het)" ~ "SDAxHET", 
    cross == "CdKO(hom)xSDA" ~ "HOMxSDA", 
    cross == "SDAxCdKO(hom)" ~ "SDAxHOM", 
    .default = as.character(cross))) %>% 
    drop_na(cross)

# remove wildtype crosses
fertility_data <- fertility_data %>% 
  filter(!cross %in% c("WTxSDA", "SDAxWT", "SDAxSDA"))

# Fecundity ====
egg_model <- glmmTMB(eggs ~ cross + line + (1|line:rep/plate_well), 
                     family =nbinom2(), 
                     ziformula=~ cross + line,
                     data = fertility_data)

emmeans::emmeans(egg_model, ~ cross + line,
                 component = "response")

#Fertility ====

# Almost all 2360 crosses produced no eggs so hatching model fails to converge
fertility_data_no2360 <-  fertility_data %>% filter(line != "2360B5")

larvae_model <- glmmTMB(cbind(larvae,(eggs-larvae)) ~ line * cross + (1|line:rep/plate_well),
                        family = binomial,
                        data = fertility_data_no2360)


# Create a combined cross and line category to analyse the 2360 group
fertility_data_combined <- fertility_data %>% 
  filter(line %in% c("1759", "QA383P") | (line == "2360B5" & cross == "SDAxHET")) %>% 
  unite("line_cross", line,cross)

larvae_model2 <- glmmTMB(cbind(larvae,(eggs-larvae)) ~ line_cross + (1|line_cross:rep/plate_well),
                        family = binomial,
                        data = fertility_data_combined)


