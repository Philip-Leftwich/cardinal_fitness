source(here::here("scripts", "packages.R"))
source(here::here("scripts", "functions", "survival_helpers.R"))

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

line_colours <- c(
  "2360B5" = "#8BABD3",
  "1759" = "#BB8BD3",
  "D251" = "#FFA040",
  "QA383P" = "#FFE699",
  "SDA-500" = "lightgray"
)

# Shared theme for transhet panels
transhet_theme <- list(
  scale_y_continuous(limits = c(0, 1)),
  labs(
    x = "Time (Hours)",
    y = "Survival Probability"
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

#' Time ratios for chosen level-versus-reference pairs of one factor
#'
#' Computes linear contrasts between levels of a single factor in a fitted
#' `flexsurvreg` accelerated failure time model (e.g. Weibull), using the
#' model coefficients and their variance-covariance matrix. This keeps the
#' shared shape parameter of the original model while allowing any level to
#' serve as the comparator (e.g. each line's own wildtype siblings), rather
#' than only the factor's reference level.
#'
#' @param model A fitted `flexsurvreg` AFT model.
#' @param factor_name Character. Name of the factor in the model formula
#'   (e.g. `"line_treatment"`).
#' @param pairs A data frame with columns `level` (the group of interest) and
#'   `ref` (its comparator), plus optionally `family` (contrasts sharing a
#'   family are adjusted together; defaults to `ref`).
#' @param adjust Character. Multiplicity adjustment passed to
#'   [stats::p.adjust()] within each family. Default `"holm"`.
#' @return `pairs` with added columns `tr`, `tr_lcl`, `tr_ucl` (time ratio
#'   and 95% CI), `z`, `p_raw` and `p_adj`. Pairs involving a level that is
#'   not in the fitted model (e.g. a group with no deaths that was fitted out
#'   to avoid complete separation) are returned with `NA` estimates.
contrast_levels <- function(model, factor_name, pairs, adjust = "holm") {
  b <- coef(model)
  V <- vcov(model)
  model_levels <- model$covdata$xlev[[factor_name]]
  coef_name <- function(lvl) paste0(factor_name, lvl)

  lin_comb <- function(level, ref) {
    if (!all(c(level, ref) %in% model_levels)) {
      return(c(est = NA_real_, se = NA_real_))
    }
    a <- setNames(numeric(length(b)), names(b))
    # The model's reference level has no coefficient (its effect is zero)
    if (coef_name(level) %in% names(b)) a[coef_name(level)] <- 1
    if (coef_name(ref) %in% names(b)) a[coef_name(ref)] <- a[coef_name(ref)] - 1
    c(est = sum(a * b), se = sqrt(drop(t(a) %*% V %*% a)))
  }

  if (!"family" %in% names(pairs)) pairs$family <- pairs$ref

  pairs |>
    mutate(
      lc = map2(level, ref, lin_comb),
      log_tr = map_dbl(lc, "est"),
      se = map_dbl(lc, "se"),
      tr = exp(log_tr),
      tr_lcl = exp(log_tr - qnorm(0.975) * se),
      tr_ucl = exp(log_tr + qnorm(0.975) * se),
      z = log_tr / se,
      p_raw = 2 * pnorm(-abs(z))
    ) |>
    mutate(p_adj = p.adjust(p_raw, method = adjust), .by = family) |>
    select(-lc)
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
#'   stripped of its variable-name prefix and a time-zero row (estimate = 1)
#'   prepended per stratum, so step lines start at 1 and groups with few or no
#'   deaths are still drawn.
tidy_km <- function(data, formula) {
  km_fit <- survfit(formula, data = data)
  strata_prefix <- paste0(all.vars(formula[[3]]), "=")
  tidy_df <- broom::tidy(km_fit) |>
    mutate(
      strata = str_remove_all(strata, paste(strata_prefix, collapse = "|"))
    )
  time_zero <- tidy_df |>
    distinct(strata) |>
    mutate(time = 0, estimate = 1)
  bind_rows(time_zero, tidy_df) |>
    arrange(strata, time)
}

#' Build a single trans-het survival panel
#'
#' Filters Weibull and Kaplan-Meier data to a single line and produces a
#' ggplot panel with ribbon CIs, model lines, and KM step overlays. Colour
#' identifies genotype; linetype identifies the estimate (solid = Weibull
#' model, dashed = Kaplan-Meier), matching the other survival figures.
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
#' @param legend_order Integer. Order of this panel's colour legend when
#'   panels are combined with collected guides; the shared estimate legend is
#'   always placed after every panel's colour legend. Default `1`.
#' @return A `ggplot` object.
plot_transhet_panel <- function(
  surv_df,
  km_df,
  line_val,
  legend_title,
  legend_labels,
  plot_title,
  colours,
  show_y_label = TRUE,
  legend_order = 1
) {
  p <- ggplot(
    filter(surv_df, line == line_val),
    aes(
      x = time,
      y = est,
      colour = treatment,
      fill = treatment
    )
  ) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = lt_model), linewidth = 0.8) +
    geom_step(
      data = filter(km_df, line == line_val),
      aes(x = time, y = estimate, colour = treatment, linetype = lt_km),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = colours,
      name = legend_title,
      labels = legend_labels
    ) +
    scale_fill_manual(
      values = colours,
      name = legend_title,
      labels = legend_labels
    ) +
    estimate_linetype_layers(colour_order = legend_order, linetype_order = 99) +
    ggtitle(plot_title) +
    transhet_theme +
    theme(legend.position = "bottom")

  if (!show_y_label) {
    p <- p + labs(y = NULL)
  }
  p
}
