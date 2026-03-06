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
#models <- lapply(dists, fit_model)
#names(models) <- dists

# Remove failed models
#models <- models[!sapply(models, is.null)]

# Extract model comparison metrics
#model_comparison <- bind_rows(lapply(models, function(mod) {
#  data.frame(
#    AIC = AIC(mod),
#    BIC = BIC(mod),
#    LogLik = logLik(mod)
#  )
#}))
#model_comparison$names <- dists
# Print model comparison
#print(model_comparison)

# model <- flexsurvreg(Surv(`Hours`, event) ~  treatment + line,
#                       data = data_long, dist = "weibull")

model <- flexsurvreg(Surv(`Hours`, event) ~  line_treatment,
                     data = data_long, dist = "weibull")

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# Hazard ratios (exponentiated coefficients)
model_summary

# tab_model(model)

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

# Kaplan-Meier estimates for overlay
km_fit <- survfit(Surv(Hours, event) ~ line_treatment, data = data_long)

km_df <- broom::tidy(km_fit) |>
  mutate(strata = str_remove(strata, "line_treatment="))  |> 
  separate(strata, into = c("line", "treatment"), sep = "_",
           extra = "merge")

# Visualisation: Weibull model ribbons + KM step lines
survival_after_blood <- ggplot(surv_df, aes(x = time, y = est,
                    fill = factor(treatment))) +
  # Model 95% CI ribbon
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  # Weibull model predicted line
  geom_line(aes(colour = factor(treatment))) +
  # KM empirical step line
  geom_step(
    data = km_df,
    aes(x = time, y = estimate, colour = factor(treatment)),
    linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time (Hours)",
       y = "Survival Probability",
       colour = "Treatment",
       fill = "Treatment",
       caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier") +
  scale_colour_brewer(palette = "Dark2", name = "Treatment") +
  scale_fill_brewer(palette = "Dark2", name = "Treatment") +
  facet_wrap(~line) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )


## TransHet survival plot

 data_transhet <- read_excel("clean_data/TransHetSurvivalData.xlsx") %>% 
  pivot_longer(cols = `2360Het`:`D251:2360WT`, 
               names_to = "line_treatment", 
               values_to = "event") %>% 
  drop_na(event) 

data_transhet <- data_transhet |>
  filter(Hours > 0)

model_transhet <- flexsurvreg(Surv(`Hours`, event) ~  line_treatment,
                     data = data_transhet, dist = "weibull")

model_summary_transhet <- tidy(model_transhet, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# Hazard ratios (exponentiated coefficients)
model_summary_transhet



# Generate predicted survival curves
pred_data_transhet <- expand_grid(
  line_treatment = unique(data_transhet$line_treatment)
)

# Calculate survival predictions
surv_pred_trans <- summary(model_transhet, newdata = pred_data_transhet, 
                     type = "survival", 
                     t = seq(0, max(data_long$Hours), length.out = 108),
                     tidy = TRUE)

# Convert to data frame for plotting
surv_df_trans <- as_tibble(surv_pred_trans) %>% 
  mutate(
    treatment = str_extract(line_treatment, "Het|WT"),
    line      = str_remove(line_treatment, "Het|WT")
  )

# Kaplan-Meier estimates for overlay
km_fit_transhet <- survfit(Surv(Hours, event) ~ line_treatment, data = data_transhet)

km_df_transhet <- broom::tidy(km_fit_transhet) |>
  mutate(strata = str_remove(strata, "line_treatment=")) |>
  mutate(
    treatment = str_extract(strata, "Het|WT"),
    line      = str_remove(strata, "Het|WT")
  )

# Visualisation: Weibull model ribbons + KM step lines
transhet_survival <- ggplot(surv_df_trans, aes(x = time, y = est, 
                    fill = factor(treatment))) +
  # Model 95% CI ribbon
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  # Weibull model predicted line
  geom_line(aes(colour = factor(treatment))) +
  # KM empirical step line
  geom_step(
    data = km_df_transhet,
    aes(x = time, y = estimate, colour = factor(treatment)),
    linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time (Hours)", 
       y = "Survival Probability",
       colour = "Treatment",
       fill = "Treatment",
       caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier") +
  scale_colour_brewer(palette = "Dark2", name = "Treatment") +
  scale_fill_brewer(palette = "Dark2", name = "Treatment") +
  facet_wrap(~line) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )
