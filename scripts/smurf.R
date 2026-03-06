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
  filter(delivery %in% "SUCROSE")

# Model=====

smurf_model <- glm(cbind(positive, negative) ~ genotype, data = data_wide, family = binomial(link = "logit"))

drop1(smurf_model, test = "Chisq")

summary(smurf_model)

# Plot smurf predictions ====

# Get emmeans predictions on response scale
smurf_emm <- emmeans::emmeans(smurf_model, ~ genotype, type = "response") |>
  as.data.frame()

# Factor levels

genotype_order <- c("WT", "KI", "KO")

smurf_emm <- smurf_emm |>
  mutate(
    genotype = factor(genotype, levels = genotype_order)
  )

# Raw proportions
raw_smurf <- data_wide |>
  mutate(
    prop = positive / (positive + negative),
    genotype = factor(genotype, levels = genotype_order)
  )

# Plot
p_smurf <- ggplot(smurf_emm, aes(x = genotype, y = prob, colour = genotype, group = genotype)) +

   # 95% CI
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_nudge(x = .1),
    width = 0.15, linewidth = 0.8
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
    alpha = 0.3, size = 1.5, shape = 16
  ) +
   scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_brewer(palette = "Dark2", name = "Genotype") +
  labs(
    x = "Delivery",
    y = "Smurf rate",
    title = "Smurf assay: proportion of smurf-positive flies"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_smurf
