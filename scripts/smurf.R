library(tidyverse)
library(janitor)
library(survival)
library(survminer)
library(readxl)
library(patchwork)
library(gtsummary)
library(glmmTMB)
library(emmeans)
library(performance)

# Load data ====

path <- "data/smurf.xlsx"

data <- read_excel(path = path, sheet = "Smurf screening", .name_repair = janitor::make_clean_names, skip =1) %>% 
  rename("id" = "x") %>% 
  select(!contains("x")) %>% 
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
  ) %>% 
  rename_with(
    .cols = 2:25,
    .fn = ~ str_replace(., "smurf_", "")) %>% 
  pivot_longer(cols = 2:25,
               names_to = "condition",
               values_to = "values") %>% 
  separate(condition, c("smurf", "genotype", "delivery")) %>% 
  mutate(id = str_remove_all(id, "[0-9]")) %>% 
  drop_na(id)
  
  

data_wide <- data %>% 
  group_by(id, smurf, genotype, delivery) %>% 
  summarise(n = sum(values, na.rm = T)) %>% 
  pivot_wider(names_from = smurf, values_from = n) %>% 
  ungroup() %>% 
  mutate(delivery = fct_relevel(delivery, "SUCROSE")) %>% 
  mutate(genotype = fct_relevel(genotype, "WT"))

data_wide <- data_wide %>% 
  filter(delivery !="XA")

# Model=====

model <- glm(cbind(positive, negative) ~ genotype * delivery, data = data_wide, family = binomial(link = "logit"))

drop1(model, test = "Chisq")

summary(model)
