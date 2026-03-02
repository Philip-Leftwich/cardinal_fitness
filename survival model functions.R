### CompaRE AND CONTRAST PARAMETRIC SURVIVAL MODELS

# Load packages
library(survival)
library(flexsurv)
library(dplyr)

# Example simulated dataset
set.seed(123)
n <- 200
data <- data.frame(
  time = rexp(n, rate = 0.1),              # survival time
  status = rbinom(n, size = 1, prob = 0.8), # event indicator
  age = rnorm(n, mean = 65, sd = 10)        # covariate
)

# Fit Kaplan-Meier
km_fit <- survfit(Surv(time, status) ~ 1, data = data)

# List of distributions to fit
dists <- c("exponential", "weibull", "gompertz", 
           "lognormal","gengamma")

# Function to fit model
fit_model <- function(dist) {
  tryCatch(
    flexsurvreg(Surv(time, status) ~ age, data = data, dist = dist),
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


# Create time points for predictions
times <- seq(0, max(data$time), length.out = 100)

# Get predicted survival for each model
predictions <- imap(models, function(mod, name) {
  
  pred <- summary(mod, newdata = data.frame(age = 65), t = times, type = "survival")

  tibble(
    time = pred[[1]]$time,
    survival = pred[[1]]$est,
    lower = pred[[1]]$lcl,
    upper = pred[[1]]$ucl,
  ) %>%
    mutate(model = name)   # now "name" is passed in by imap
})


predictions_df <- do.call(rbind, predictions)

ggplot() +
  geom_step(aes(x = km_fit$time, y = km_fit$surv), color = "black", size = 1.2, direction = "hv") +
  geom_ribbon(data = predictions_df,
            aes(x = time, ymin = lower, ymax = upper, fill = model), alpha = 0.2,
            linewidth = 1) +
  labs(title = "Kaplan-Meier vs Parametric Survival Models",
       x = "Time", y = "Survival Probability") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")



# Optionally: plot survival curves for visual comparison
plot(models[[1]], main = paste("Fitted Survival: ", models[[1]]$dist), col = 1)
for (i in 2:length(models)) {
  lines(models[[i]], col = i)
}
