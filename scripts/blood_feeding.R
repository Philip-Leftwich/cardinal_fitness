library(tidyverse)
library(survival)
library(survminer)
library(readxl)
library(patchwork)
library(gtsummary)
library(flexsurv)
library(gghighlight)

## Survival after BF====
path <- "clean_data/SurvivalAfterBFData.xlsx"

data <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path, sheet = .x), .id = "line")


data_long <- data %>% 
  pivot_longer(cols = `WT`:`Hom`, 
               names_to = "treatment", 
               values_to = "event") %>% 
  drop_na(event) %>% 
#  filter(line %in% c("1759", "2360", "2072")) %>% 
  unite("line_treatment", c(line,treatment),
        remove = FALSE)



# Function to fit model
fit_model <- function(dist) {
  tryCatch(
    flexsurvreg(Surv(`Hours`, event) ~ treatment + line, data = data_long, dist = dist),
    error = function(e) NULL  # return NULL if model fails
  )
}


# Fit all models
models <- lapply(dists, fit_model)
names(models) <- dists

# Remove failed models
models <- models[!sapply(models, is.null)]

# Extract model comparison metrics
model_comparison <- bind_rows(lapply(models, function(mod) {
  data.frame(
    AIC = AIC(mod),
    BIC = BIC(mod),
    LogLik = logLik(mod)
  )
}))
model_comparison$names <- dists
# Print model comparison
print(model_comparison)

# model <- flexsurvreg(Surv(`Hours`, event) ~  treatment + line,
#                       data = data_long, dist = "weibull")

model <- flexsurvreg(Surv(`Hours`, event) ~  line_treatment,
                     data = data_long, dist = "weibull")

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# Hazard ratios (exponentiated coefficients)
model_summary

# Generate predicted survival curves
pred_data <- expand_grid(
  line_treatment = unique(data_long$line_treatment)
)

# Calculate survival predictions
surv_pred <- summary(model, newdata = pred_data, 
                     type = "survival", 
                     t = seq(0, max(data_long$Hours), length.out = 108),
                     tidy = TRUE)

# Convert to data frame for plotting
surv_df <- as_tibble(surv_pred) %>% 
  separate(line_treatment,
           c("line", "treatment"))

# Visualisation
ggplot(surv_df, aes(x = time, y = est, 
                    fill = factor(treatment),
                    linetype = line)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2)+
  geom_line(aes(colour = factor(treatment))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time (Hours)", 
       y = "Survival Probability",
       colour = "Treatment",
       fill = "Treatment",
       linetype = "Genotype") +
  facet_wrap(~line)+
  theme_minimal()


## TransHet survival plot

data_transhet <- data %>% 
  pivot_longer(cols = `WT`:`Hom`, 
               names_to = "treatment", 
               values_to = "event") %>% 
  drop_na(event) %>% 
  filter(line %in% c("2360.2", "D251.2360")) %>% 
  unite("line_treatment", c(line,treatment),
        remove = FALSE)


model <- flexsurvreg(Surv(`Hours`, event) ~  line_treatment,
                     data = data_transhet, dist = "weibull")

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# Hazard ratios (exponentiated coefficients)
model_summary

# Generate predicted survival curves
pred_data <- expand_grid(
  line_treatment = unique(data_transhet$line_treatment)
)

# Calculate survival predictions
surv_pred <- summary(model, newdata = pred_data, 
                     type = "survival", 
                     t = seq(0, max(data_long$Hours), length.out = 108),
                     tidy = TRUE)

# Convert to data frame for plotting
surv_df <- as_tibble(surv_pred) %>% 
  separate(line_treatment,
           c("line", "treatment"), sep ="_")

# Visualisation
ggplot(surv_df, aes(x = time, y = est, 
                    fill = factor(treatment),
                    linetype = line)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2)+
  geom_line(aes(colour = factor(treatment))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time (Hours)", 
       y = "Survival Probability",
       colour = "Treatment",
       fill = "Treatment",
       linetype = "Genotype") +
  facet_wrap(~line)+
  theme_minimal()
