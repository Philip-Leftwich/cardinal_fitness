source(here::here("scripts", "packages.R"))

# ── Shared lookups ─────────────────────────────────────────────────────────────

line_labels <- c(
  "1759" = expression(italic(cd)^g384),
  "2360B5" = expression(italic(cd)^g225),
  "QA383P" = expression(
    italic(cd)^{
      "384R"
    }
  )
)

line_colours <- c(
  "2360B5" = "#8BABD3",
  "1759" = "#BB8BD3",
  "D251" = "#FFA040",
  "QA383P" = "#FFE699",
  "SDA-500" = "darkgray"
)

line_shape <- c("1759" = 16, "2360B5" = 17, "QA383P" = 15)

# Shared plot theme
fertility_theme <- list(
  scale_colour_manual(
    values = line_colours,
    name = "Line",
    labels = line_labels
  ),
  scale_shape_manual(values = line_shape, name = "Line", labels = line_labels),
  labs(x = "Cross"),
  theme_bw(base_size = 12),
  theme(
    axis.text.x = element_text(hjust = 0.5),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )
)

# ── Helper functions ───────────────────────────────────────────────────────────

#' Extract emmeans predictions on the response scale
#'
#' Computes marginal means from a fitted model over `line` and `cross`,
#' back-transforms to the response scale, and re-levels `cross`.
#'
#' @param model A fitted model accepted by [emmeans::emmeans()] (e.g.
#'   a `glmmTMB` object).
#' @param cross_levels Character vector. Factor levels for `cross` in the
#'   desired display order.
#' @return A tibble of emmeans predictions with `NA` rows dropped and `cross`
#'   re-levelled to `cross_levels`.
get_emm <- function(model, cross_levels) {
  emmeans::emmeans(model, ~ line + cross, type = "response") |>
    as_tibble() |>
    mutate(cross = factor(cross, levels = cross_levels)) |>
    drop_na()
}

get_emm_contrasts <- function(model) {
  emm1 <- emmeans::emmeans(model, pairwise ~ line | cross)
  contrasts1 <- as_tibble(emm1$contrasts) |>
    drop_na() |>
    mutate(
      ratio = exp(estimate),
      ratio.LCL = exp(estimate - 1.96 * SE),
      ratio.UCL = exp(estimate + 1.96 * SE)
    )

  emm2 <- emmeans::emmeans(model, pairwise ~ cross | line)
  contrasts2 <- as_tibble(emm2$contrasts) |>
    drop_na() |>
    mutate(
      ratio = exp(estimate),
      ratio.LCL = exp(estimate - 1.96 * SE),
      ratio.UCL = exp(estimate + 1.96 * SE)
    )

  list(
    line_by_cross = contrasts1,
    cross_by_line = contrasts2
  )
}


#' Plot model predictions and raw data for a fertility outcome
#'
#' Builds a ggplot showing emmeans point estimates with asymptotic 95\% CIs,
#' dodged by `line`, overlaid with jittered raw observations.
#'
#' @param emm_data A tibble of emmeans predictions from [get_emm()] with
#'   columns `cross`, `response`, `line`, `asymp.LCL`, `asymp.UCL`.
#' @param raw_data A data frame of raw observations to overlay as jittered
#'   points.
#' @param raw_y <tidy-select> The unquoted column in `raw_data` to use as the
#'   y aesthetic for raw points.
#' @param y_label Character. Y-axis label.
#' @param plot_title Character. Plot title.
#' @param y_scale A ggplot2 scale object to apply to the y axis (e.g.
#'   `scale_y_continuous(...)`), or `NULL` to use the default. Default `NULL`.
#' @return A `ggplot` object.
plot_fertility <- function(
  emm_data,
  raw_data,
  raw_y,
  y_label,
  y_scale = NULL
) {
  p <- ggplot(
    emm_data,
    aes(x = cross, y = response, colour = line, group = line, shape = line)
  ) +
    geom_point(
      data = raw_data,
      aes(y = {{ raw_y }}, colour = line),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
      alpha = 0.3,
      size = 1.5,
      shape = 16
    ) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    geom_errorbar(
      aes(ymin = asymp.LCL, ymax = asymp.UCL),
      position = position_dodge(width = 0.6),
      width = 0.25,
      linewidth = 0.8
    ) +
    labs(y = y_label) +
    fertility_theme

  if (!is.null(y_scale)) {
    p <- p + y_scale
  }
  p
}
