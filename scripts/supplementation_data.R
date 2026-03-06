library(tidyverse)
library(survival)
library(survminer)
library(readxl)
library(patchwork)
library(sjPlot)
library(flexsurv)

# Load data ====

path <- "clean_data/XASupplementationData.xlsx"

data <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path, sheet = .x), .id = "experiment") %>% 
  filter(experiment != "Additional assays")

# Shape data ====
data_long <- data %>% 
  pivot_longer(cols = `WT (0mM)`:`Het (24mM)`, 
               names_to = "treatment", 
               values_to = "event") %>% 
  drop_na(event)

data_long



# Format data for interaction models

data_long2 <- data_long %>%
  separate(treatment, into = c("genotype", "dosage"), sep = "\\(") %>% 
  mutate(genotype = str_trim(genotype)) %>% 
  mutate(dosage = str_remove(dosage, "\\)")) %>% 
  mutate(experiment = str_remove(experiment, "XA ")) %>% 
  separate(dosage, c("dosage", "c-")) %>% 
  mutate("c_minus" = factor(if_else(`c-` == "C", "-", "+", missing = "+"))) %>% 
  mutate(dosage = factor(str_remove(dosage, "mM"),
                         levels = c("0","6","24"))) %>% 
  rename("dosage_mM" = dosage)

# coxph
# models <- data_long2 %>% 
#   group_by(supplementation) %>% 
#   nest() %>% 
#   mutate(model = map(data, ~coxph(Surv(`Hours`, event) ~ source + dosage + genotype, data = .) %>% 
#                        tbl_regression(exp = TRUE)))
# 
# 
# models <- data_long2 %>% 
#   group_by(supplementation) %>% 
#   nest() %>% 
#   mutate(data = map(data, ~{
#     .x %>% mutate(dosage_mM = fct_drop(dosage_mM))})) %>% 
#   mutate(model = map(data, ~survfit(Surv(`Hours`, event) ~ source + dosage_mM + genotype, data = .))) %>% 
#   mutate(plot = map2(.x = model, .y = data, ~ggsurvplot(.x, data = .y)))
# 
# walk(models$plot, ~print(.x))


## Parametric  models

dists <- c("exponential", "weibull", "gompertz", 
           "lognormal","gengamma")

# Function to fit model
fit_model <- function(data, dist) {
  tryCatch(
    flexsurvreg(Surv(`Hours`, event) ~ experiment + dosage_mM + genotype + experiment:genotype, data = data, dist = dist),
    error = function(e) NULL  # return NULL if model fails
  )
}


test <- data_long2 %>% 
  nest() %>% 
  mutate(data = map(data, ~{
    .x %>% mutate(dosage_mM = fct_drop(dosage_mM))})) %>% 
  expand_grid(distribution = dists) %>% 
  mutate(model = map2(data, distribution, ~fit_model(.x, .y))) 
                       
nested_models <- test %>%
  mutate(
    aic = map_dbl(
      model,
      ~ if (!is.null(.x)) AIC(.x) else NA_real_
    )
  ) %>%  
  arrange(aic)


# weibull model====

model <- flexsurvreg(Surv(`Hours`, event) ~ experiment + 
                       dosage_mM + 
                       genotype + 
                       `c_minus`+
                       experiment:dosage_mM +
                       experiment:genotype + 
                       experiment:`c_minus` +
                       dosage_mM:genotype +
                       genotype:`c_minus`,
                     data = data_long2, dist = "weibull")

model2 <- flexsurvreg(Surv(`Hours`, event) ~ experiment + 
                       dosage_mM + 
                       genotype + 
                       `c_minus`+
                       experiment:dosage_mM +
                       experiment:genotype + 
                       experiment:`c_minus`,
                     data = data_long2, dist = "weibull")

model3 <- flexsurvreg(Surv(`Hours`, event) ~ experiment + 
                        dosage_mM + 
                        genotype + 
                        `c_minus`+
                        experiment:dosage_mM +
                        experiment:`c_minus`,
                      data = data_long2, dist = "weibull")

model4 <- flexsurvreg(Surv(`Hours`, event) ~ experiment + 
                        dosage_mM + 
                        genotype+ 
                        c_minus+
                        experiment:c_minus,
                      data = data_long2, dist = "weibull")

# Extract model coefficients with confidence intervals
model_summary <- tidy(model4, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# Hazard ratios (exponentiated coefficients)
model_summary

model4

# Generate predicted survival curves
pred_data <- expand_grid(
  experiment = unique(data_long2$experiment),
  dosage_mM = unique(factor(c(0,6,24))),
  genotype = unique(data_long2$genotype),
  c_minus = unique(data_long2$c_minus)
)

# Calculate survival predictions
surv_pred <- summary(model4, newdata = pred_data, 
                     type = "survival", 
                     t = seq(0, max(data_long2$Hours), length.out = 108),
                     tidy = TRUE)

# Convert to data frame for plotting
surv_df <- as_tibble(surv_pred)

# Visualisation
XA_suppl <- ggplot(surv_df, aes(x = time, y = est, 
                    fill = factor(dosage_mM),
                    linetype = genotype)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2)+
  geom_line(aes(colour = factor(dosage_mM))) +
  facet_wrap(c_minus~experiment) +
 scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time (Hours)", 
       y = "Survival Probability",
       colour = "Treatment",
       fill = "Treatment") +
  scale_colour_manual(values = c("#FBEAFF", "#B39CD0", "#845EC2"), name = "Dosage")+
  scale_fill_manual(values = c("#FBEAFF", "#B39CD0", "#845EC2"), name = "Dosage")+
   theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )
# Model fit statistics
glance(model4)

# ### Predictions====
# 
# test_2 <- test %>% filter(distribution == "weibull")
# 
# times <- seq(0, max(data_long2$Hours, na.rm = TRUE), length.out = 200)
# 
# 
# predictions <- data_long2 %>% 
#   select(experiment, dosage_mM, genotype) %>% 
#   distinct() %>% 
#   nest()
# 
# 
# 
# test_3 <- test_2 %>%
#   mutate(predictions = map2(model, predictions$data, function(mod, nested_df) {
#     
#     nested_df %>%
#       pmap_dfr(function(experiment,dosage_mM, genotype) {
#         newdata <- data.frame(
#           experiment = experiment,
#           source = source,
#           dosage_mM = dosage_mM,
#           genotype = genotype
#         )
#         
#         # Try safely getting the prediction summary
#         pred <- tryCatch(
#           summary(mod, newdata = newdata, t = times, type = "survival"),
#           error = function(e) NULL
#         )
#         
#         # Check that pred is valid and has expected structure
#         if (is.list(pred) && !is.null(pred[[1]]$time)) {
#           pred_df <- pred[[1]]
#           tibble(
#             time = pred_df$time,
#             survival = pred_df$est,
#             lower = pred_df$lcl,
#             upper = pred_df$ucl,
#             dosage_mM = dosage_mM,
#             genotype = genotype,
#             experiment = experiment
#           )
#         } else {
#           tibble(
#             time = NA,
#             survival = NA,
#             lower = NA,
#             upper = NA,
#             dosage_mM = dosage_mM,
#             genotype=genotype,
#             experiment = experiment
#           )
#         }
#       })
#   }))
# 
# # Fit Kaplan-Meier
# km_fit <- survfit(Surv(`Hours`, event) ~ source + dosage + genotype, data = test[[10,2]][[1]])
# 
# km_fit <- broom::tidy(km_fit) %>% rename(group = strata)
# 
# predictions_df <- predictions_df %>%  rename(group = combo_label)
# 
# ggplot() +
#   geom_step( data = km_fit, aes(x = time, y = estimate, colour = group), size = 1.2, direction = "hv") +
#   geom_ribbon(data = predictions_df,
#               aes(x = time, ymin = lower, ymax = upper, fill = group), alpha = 0.2,
#               linewidth = 1) +
#   labs(title = "Kaplan-Meier vs Parametric Survival Models",
#        x = "Time", y = "Survival Probability") +
#   theme_minimal() +
#   guides(fill = "none")+
#   scale_fill_brewer(palette = "Set1")
# 
# # Fit weibull
#  
#     
#  
#     plots <- map2(test_3$predictions, test_3$supplementation, ~ ggplot(data = .x,
#     aes(x = time, y = survival)) +
#     geom_smooth(aes(linetype = genotype, 
#                     colour = dosage_mM,
#                     group = interaction(genotype, dosage_mM)), 
#               direction = "hv") +
#     geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = dosage_mM, group = interaction(genotype,dosage_mM)), alpha = 0.2,
#                 linewidth = 1) +
#     # gghighlight()+
#     labs(x = "Time", y = "Survival Probability") +
#     theme_minimal() +
#     #  scale_fill_brewer(palette = "Set1")+
#     #  scale_colour_brewer(palette = "Set1")+
#     guides(colour  = "none")+
#   #  scale_linetype_manual(values=c("twodash", "dotted", "solid"))+
#     scale_y_continuous(labels = scales::percent, limits = c(0,1))+
#     theme(legend.key.size = unit(0.5, "in"))+
#     facet_wrap(~source)+
#     scale_fill_discrete_sequential(palette = "Rocket")+
#     scale_color_discrete_sequential(palette = "Rocket")+
#     ggtitle(glue::glue("{.y} supplementation"))
#     )
# 
# # Model tables
#     summarize_flexsurv(test_3$model[[1]]) %>% 
#       gt() %>% 
#       tab_header(
#         title = md(glue::glue("{test_3$supplementation[[1]]}")))
#     
#         
#     summarize_flexsurv(test_3$model[[2]]) %>% 
#       gt() %>% 
#       tab_header(
#         title = md(glue::glue("{test_3$supplementation[[2]]}")))
# 
#     summarize_flexsurv(test_3$model[[3]]) %>% 
#       gt() %>% 
#       tab_header(
#         title = md(glue::glue("{test_3$supplementation[[3]]}")))
# 
# # Split data by supplementation into a named list
# data_list <- data_long2 %>% split(.$supplementation)
# 
# # Initialize empty lists to store models and plots
# model_list <- list()
# plot_list <- list()
# 
# # Loop over each supplementation group
# for (group_name in names(data_list)) {
#   
#   # Extract the data for this group
#   df <- data_list[[group_name]]
#   
#   # Fit the model
#   fit <- flexsurvreg(Surv(Hours, event) ~ source + dosage + genotype, data = df, dist = "gengamma")
#   
#   # Store the model in the list using group name as key
#   model_list[[group_name]] <- fit
#   
#   # Generate the plot and store it
#   plot_list[[group_name]] <- ggflexsurvplot(fit = fit, data = df)
# }
# 
# 
# 
# 
# 
# 
# ## Trp
# 
# trp_model <- flexsurvreg(Surv(`Hours`, event) ~ source + dosage + genotype + source:dosage + dosage:genotype + source:genotype, data = Trp, dist = "gengamma")
# 
# # best AIC
# trp_model <- flexsurvreg(Surv(`Hours`, event) ~ source + dosage + genotype  + source:genotype, data = Trp, dist = "gengamma")
# 
# cox_model <- coxph(Surv(`Hours`, event) ~ source + dosage + genotype  + source:genotype, data = Trp)
# 
# zph <- cox.zph(cox_model)
# 
# # Then plot with ggcoxzph()
# ggcoxzph(zph)
# 
# 
# 
# # Trp <- Trp %>% mutate(dose_numeric = as.numeric(str_remove(dosage, "mM")))
# # OK? 
# 
# #XA <- XA %>% separate(dosage, into = c("dose", "add"), sep = " ") %>% mutate(add = if_else(is.na(add), "C+",add)) %>% mutate(dose_numeric = as.numeric(str_remove(dose, "mM")))
# 
# cox_model_2 <- coxph(Surv(`Hours`, event) ~ strata(source) + dosage * genotype, data = Trp_2)
# 
# 
# cox_model_3 <- coxph(Surv(`Hours`, event) ~ strata(source) + dose_numeric + genotype, data = Trp_2)
# 
# cox_model_4 <- coxph(Surv(`Hours`, event) ~ source + dose_numeric + genotype, data = Gly)
# 
# cox_model_5 <- coxph(Surv(`Hours`, event) ~ strata(source) + dose_numeric + add + genotype, data = XA)


