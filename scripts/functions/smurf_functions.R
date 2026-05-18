source(here::here("scripts", "packages.R"))

# ── Helper functions ───────────────────────────────────────────────────────────

#' Fit a Weibull survival model and return predicted survival curves
#'
#' Fits a parametric Weibull model via [flexsurv::flexsurvreg()] and generates
#' a tidy tibble of predicted survival probabilities over a time grid,
#' automatically building the prediction grid from unique covariate values.
#'
#' @param data A data frame containing the survival data.
#' @param formula A formula passed to [flexsurv::flexsurvreg()], typically
#'   `Surv(time, event) ~ covariates`.
#' @param t_max Numeric. Upper bound of the prediction time sequence.
#' @param n_t Integer. Number of evenly-spaced time points to predict at.
#'   Default `108`.
#' @return A named list with two elements:
#'   \describe{
#'     \item{model}{The fitted `flexsurvreg` object.}
#'     \item{surv_df}{A tibble of predicted survival probabilities with columns
#'       `time`, `est`, `lcl`, `ucl`, and one column per covariate.}
#'   }
fit_weibull_surv <- function(data, formula, t_max, n_t = 108) {
  model <- flexsurvreg(formula, data = data, dist = "weibull")
  model$call$data <- substitute(data)
  model$call$formula <- substitute(formula)
  pred_vars <- all.vars(formula)[-c(1, 2)] # predictor names from formula
  pred_data <- map(pred_vars, ~ unique(data[[.x]])) |>
    setNames(pred_vars) |>
    do.call(what = expand_grid, args = _)
  surv_pred <- summary(
    model,
    newdata = pred_data,
    type = "survival",
    t = seq(0, t_max, length.out = n_t),
    tidy = TRUE
  )
  list(model = model, surv_df = as_tibble(surv_pred))
}


#' Fit and tidy a Kaplan-Meier survival curve
#'
#' Fits a Kaplan-Meier estimator, strips the variable-name prefix from strata
#' labels (e.g. `"group=A"` -> `"A"`), and prepends a time-zero row
#' (estimate = 1) for each stratum so step lines start at 1.
#'
#' @param data A data frame containing the survival data.
#' @param formula A formula passed to [survival::survfit()], typically
#'   `Surv(time, event) ~ strata_var`.
#' @return A tibble from [broom::tidy.survfit()] with cleaned `strata` labels
#'   and an added time-zero row per stratum.
tidy_km <- function(data, formula) {
  km_fit <- survfit(formula, data = data)
  strata_prefix <- paste0(all.vars(formula[[3]]), "=")
  tidy_df <- broom::tidy(km_fit) |>
    mutate(across(
      strata,
      ~ str_remove_all(.x, paste(strata_prefix, collapse = "|"))
    ))
  # Prepend time = 0, estimate = 1 for each stratum so step lines start at 1
  time_zero <- tidy_df |>
    distinct(strata) |>
    mutate(time = 0, estimate = 1)
  bind_rows(time_zero, tidy_df) |>
    arrange(strata, time)
}

#' Build a Weibull + Kaplan-Meier survival plot
#'
#' Combines a shaded ribbon CI and model line from Weibull predictions with a
#' dashed KM step overlay.
#'
#' @param surv_df A tibble of Weibull predicted survival from
#'   [fit_weibull_surv()] with columns `time`, `est`, `lcl`, `ucl`, and a
#'   grouping column named by `colour_var`.
#' @param km_df A tibble of tidy KM estimates from [tidy_km()] with columns
#'   `time`, `estimate`, and a strata column named by `colour_var`.
#' @param colour_var Character. Name of the column used for colour/fill
#'   grouping in both `surv_df` and `km_df`.
#' @param colour_scale A ggplot2 colour scale (e.g.
#'   `scale_colour_brewer(...)`).
#' @param fill_scale A ggplot2 fill scale (e.g. `scale_fill_brewer(...)`).
#' @param extra_layers A ggplot2 layer or list of layers to append (e.g.
#'   `facet_wrap(...)`), or `NULL`. Default `NULL`.
#' @return A `ggplot` object.
plot_surv <- function(
  surv_df,
  km_df,
  colour_var,
  colour_scale,
  fill_scale,
  extra_layers = NULL
) {
  p <- ggplot(surv_df, aes(x = time, y = est, fill = .data[[colour_var]])) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
    geom_line(aes(colour = .data[[colour_var]])) +
    geom_step(
      data = km_df,
      aes(x = time, y = estimate, colour = .data[[colour_var]]),
      linetype = "dashed",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 120, by = 24)) +
    colour_scale +
    fill_scale +
    labs(
      x = "Time (Hours)",
      y = "Survival Probability",
      caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank()
    )
  if (!is.null(extra_layers)) {
    p <- p + extra_layers
  }
  p
}
