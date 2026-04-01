source(here::here("scripts", "packages.R"))

# ── Shared lookups ─────────────────────────────────────────────────────────────

line_labels <- c(
  "1759" = "italic(cd)^g384",
  "2360B5" = "italic(cd)^g225",
  "2072" = "italic(cd)^g384_del",
  "2301" = "italic(cd)^{g338-384}",
  "D251" = "italic(cd)^{225*R}",
  "QA383P" = "italic(cd)^{384*R}"
)

treatment_colours <- c(
  "WT" = "#4D4D4D",
  "Het" = "#E69F00",
  "Hom" = "#56B4E9"
)

# Shared theme for transhet panels
transhet_theme <- list(
  scale_y_continuous(limits = c(0, 1)),
  labs(
    x = "Time (Hours)",
    y = "Survival Probability",
    caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier"
  ),
  theme_bw(base_size = 12),
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    panel.grid.minor = element_blank()
  )
)

# ── Helper functions ───────────────────────────────────────────────────────────

#' Fit a Weibull parametric survival model
#'
#' @param data A data frame containing the survival data.
#' @param formula A formula passed to [flexsurv::flexsurvreg()], typically
#'   of the form `Surv(time, event) ~ covariates`.
#' @return A named list with two elements:
#'   \describe{
#'     \item{model}{The fitted `flexsurvreg` object.}
#'     \item{summary}{A tibble of exponentiated coefficients (HRs) with 95\% CIs,
#'       rounded to 3 decimal places.}
#'   }
fit_weibull <- function(data, formula) {
  model <- flexsurvreg(formula, data = data, dist = "weibull")
  model$call$data <- substitute(data)
  model$call$formula <- substitute(formula)
  summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  list(model = model, summary = summary)
}

#' Predict survival probabilities from a fitted Weibull model
#'
#' @param model A fitted `flexsurvreg` object.
#' @param data The data frame used to fit `model`. Unique values of
#'   `line_treatment` are used as the prediction grid.
#' @param t_max Numeric. Upper bound of the time sequence for predictions.
#' @param n_t Integer. Number of time points to predict at. Default `108`.
#' @return A tibble of predicted survival probabilities with columns
#'   `time`, `est`, `lcl`, `ucl`, and `line_treatment`.
pred_weibull <- function(model, data, t_max, n_t = 108) {
  pred_data <- expand_grid(line_treatment = unique(data$line_treatment))
  summary(
    model,
    newdata = pred_data,
    type = "survival",
    t = seq(0, t_max, length.out = n_t),
    tidy = TRUE
  ) |>
    as_tibble()
}

#' Fit and tidy a Kaplan-Meier survival curve
#'
#' Fits a Kaplan-Meier estimator and returns a tidy tibble with strata
#' labels cleaned of their variable-name prefix (e.g. `"group=A"` -> `"A"`).
#'
#' @param data A data frame containing the survival data.
#' @param formula A formula passed to [survival::survfit()], typically
#'   of the form `Surv(time, event) ~ strata_var`.
#' @return A tibble from [broom::tidy.survfit()] with the `strata` column
#'   stripped of its variable-name prefix.
tidy_km <- function(data, formula) {
  km_fit <- survfit(formula, data = data)
  strata_prefix <- paste0(all.vars(formula[[3]]), "=")
  broom::tidy(km_fit) |>
    mutate(
      strata = str_remove_all(strata, paste(strata_prefix, collapse = "|"))
    )
}

#' Build a single trans-het survival panel
#'
#' Filters Weibull and Kaplan-Meier data to a single line and produces a
#' ggplot panel with ribbon CIs, model lines, and KM step overlays.
#'
#' @param surv_df A tibble of Weibull predicted survival (from [pred_weibull()]).
#'   Must contain columns `line`, `time`, `est`, `lcl`, `ucl`, `treatment`.
#' @param km_df A tibble of tidy KM estimates (from [tidy_km()]).
#'   Must contain columns `line`, `time`, `estimate`, `treatment`.
#' @param line_val Character. The value of `line` to filter to for this panel.
#' @param legend_title Character. HTML/markdown string used as the legend title.
#' @param legend_labels Named character vector mapping treatment levels to
#'   display labels (supports HTML via `ggtext`).
#' @param plot_title An expression or character string for the panel title.
#' @param show_y_label Logical. If `FALSE`, the y-axis label is suppressed
#'   (useful for right-hand panels in a patchwork). Default `TRUE`.
#' @return A `ggplot` object.
plot_transhet_panel <- function(
  surv_df,
  km_df,
  line_val,
  legend_title,
  legend_labels,
  plot_title,
  show_y_label = TRUE
) {
  p <- ggplot(
    filter(surv_df, line == line_val),
    aes(
      x = time,
      y = est,
      colour = treatment,
      linetype = treatment,
      fill = treatment
    )
  ) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_step(
      data = filter(km_df, line == line_val),
      aes(x = time, y = estimate, colour = treatment, linetype = treatment),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = c("Het" = "#E69F00", "WT" = "#4D4D4D"),
      name = legend_title,
      labels = legend_labels
    ) +
    scale_fill_manual(
      values = c("Het" = "#E69F00", "WT" = "#4D4D4D"),
      name = legend_title,
      labels = legend_labels
    ) +
    scale_linetype_manual(
      values = c("Het" = "solid", "WT" = "dashed"),
      name = legend_title,
      labels = legend_labels
    ) +
    ggtitle(plot_title) +
    transhet_theme

  if (!show_y_label) {
    p <- p + labs(y = NULL)
  }
  p
}
