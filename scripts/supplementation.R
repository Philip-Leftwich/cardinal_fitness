library(tidyverse)
library(survival)
library(survminer)
library(readxl)
library(patchwork)
library(gtsummary)

# Load data ====

path <- "data/supplementation-data.xlsx"

data <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path, sheet = .x), .id = "experiment") %>% 
  filter(experiment != "Additional assays")

# Shape data ====
data_long <- data %>% 
  pivot_longer(cols = `WT (0mM)`:`Het (50mM)`, 
               names_to = "treatment", 
               values_to = "event") %>% 
  drop_na(event)

data_long


# Recreate original plots ====

models <- data_long %>% 
  group_by(experiment) %>% 
  nest() %>% 
  mutate(model = map(data, ~survfit(Surv(`Hours`, event) ~ treatment, data = .))) %>% 
  mutate(plot = map2(.x = model, .y = data, ~ggsurvplot(.x, data = .y)))

walk(models$plot, ~print(.x))

# Format data for interaction models

data_long2 <- data_long %>%
  separate(treatment, into = c("genotype", "dosage"), sep = "\\(") %>% 
  mutate(genotype = str_trim(genotype)) %>% 
  mutate(dosage = str_remove(dosage, "\\)")) %>% 
  separate(experiment, into = c("supplementation", "source"), remove = FALSE, extra = "merge")


models <- data_long2 %>% 
  group_by(supplementation) %>% 
  nest() %>% 
  mutate(model = map(data, ~coxph(Surv(`Hours`, event) ~ source + dosage + genotype, data = .) %>% 
                       tbl_regression(exp = TRUE)))



models <- data_long2 %>% 
  group_by(supplementation) %>% 
  nest() %>% 
  mutate(model = map(data, ~survfit(Surv(`Hours`, event) ~ source + dosage + genotype, data = .))) %>% 
  mutate(plot = map2(.x = model, .y = data, ~ggsurvplot(.x, data = .y)))

walk(models$plot, ~print(.x))
