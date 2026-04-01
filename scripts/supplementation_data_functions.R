source(here::here("scripts", "packages.R"))

# ── Helper functions ───────────────────────────────────────────────────────────

#' Safely fit a flexible parametric survival model
#'
#' Wraps [flexsurv::flexsurvreg()] in [tryCatch()] so that failed model fits
#' (e.g. non-convergence) return `NULL` rather than an error, making it safe
#' to use inside `map()`.
#'
#' @param formula A formula passed to [flexsurv::flexsurvreg()], typically
#'   `Surv(time, event) ~ covariates`.
#' @param data A data frame containing the survival data.
#' @param dist Character. Distribution name passed to `flexsurvreg()`
#'   (e.g. `"weibull"`, `"gompertz"`, `"gengamma"`).
#' @return A fitted `flexsurvreg` object, or `NULL` if the model fails.
fit_flexsurv <- function(formula, data, dist) {
  tryCatch(
    flexsurvreg(formula, data = data, dist = dist),
    error = function(e) NULL
  )
}

#' Tidy a flexsurvreg model summary
#'
#' Returns a tibble of exponentiated coefficients (time ratios / hazard ratios)
#' with 95\% CIs, rounded to 3 decimal places.
#'
#' @param model A fitted `flexsurvreg` object.
#' @return A tibble from [broom::tidy()] with numeric columns rounded to 3 dp
#'   and estimates exponentiated.
tidy_model <- function(model) {
  tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

#' Generate predicted survival curves from a flexsurvreg model
#'
#' Calls [summary.flexsurvreg()] with `tidy = TRUE` over an evenly-spaced time
#' grid and returns the result as a tibble.
#'
#' @param model A fitted `flexsurvreg` object.
#' @param newdata A data frame of covariate combinations to predict for.
#' @param t_max Numeric. Upper bound of the prediction time sequence.
#' @param n_t Integer. Number of evenly-spaced time points. Default `108`.
#' @return A tibble with columns `time`, `est`, `lcl`, `ucl`, and one column
#'   per covariate in `newdata`.
pred_surv <- function(model, newdata, t_max, n_t = 108) {
  summary(
    model,
    newdata = newdata,
    type = "survival",
    t = seq(0, t_max, length.out = n_t),
    tidy = TRUE
  ) |>
    as_tibble()
}
