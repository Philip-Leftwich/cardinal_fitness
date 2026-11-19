source(here::here("scripts", "packages.R"))

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

#' Compute pairwise emmeans contrasts for a fertility model
#'
#' Returns pairwise contrasts from a fitted model for two stratifications:
#' line comparisons within each cross type, and cross comparisons within each
#' line. Estimates are exponentiated to the ratio scale with 95% CIs.
#'
#' @param model A fitted model accepted by [emmeans::emmeans()] (e.g. a
#'   `glmmTMB` object).
#' @return A named list with two tibbles:
#'   \describe{
#'     \item{line_by_cross}{Pairwise line contrasts within each cross type.}
#'     \item{cross_by_line}{Pairwise cross contrasts within each line.}
#'   }
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


#' Contrasts of each cross against its own line's wildtype group
#'
#' Treatment-versus-control contrasts (`trt.vs.ctrl`) of every cross against
#' the `WT` cross within each line, on the ratio scale (incidence rate ratio
#' for a log-link count model; odds ratio for a logit-link binomial model).
#' p-values are Dunnett-adjusted within line; 95% CIs are unadjusted,
#' matching [get_emm_contrasts()].
#'
#' @param model A fitted model accepted by [emmeans::emmeans()] with factors
#'   `line` and `cross`.
#' @param ref Character. The control level of `cross`. Default `"WT"`.
#' @return A tibble with one row per non-WT cross within each line.
get_emm_vs_wt <- function(model, ref = "WT") {
  emm <- emmeans::emmeans(model, ~ cross | line)
  ref_idx <- which(levels(emm)$cross == ref)
  emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_idx) |>
    summary(infer = c(FALSE, TRUE)) |>
    as_tibble() |>
    # drop non-estimable cells (e.g. cd[g225] HET females, which laid no eggs)
    filter(is.finite(estimate), is.finite(SE)) |>
    mutate(
      ratio = exp(estimate),
      ratio.LCL = exp(estimate - qnorm(0.975) * SE),
      ratio.UCL = exp(estimate + qnorm(0.975) * SE)
    )
}

#' Split the combined cross label into genotype and transgenic parent
#'
#' `cross` combines the transgenic parent's sex and genotype (e.g.
#' "Female \nHET" = heterozygous transgenic females crossed to SDA-500
#' males). This separates the two so that Figure 1 can show genotype on the
#' x axis and cross direction by shape.
#'
#' @param df A data frame with a `cross` column.
#' @return `df` with factor columns `genotype` (WT, HET, HOM) and `tg_parent`
#'   (Transgenic mother, Transgenic father, Wildtype cross).
split_cross <- function(df) {
  df |>
    mutate(
      genotype = case_when(
        cross == "WT" ~ "WT",
        str_detect(cross, "HET") ~ "HET",
        str_detect(cross, "HOM") ~ "HOM"
      ) |>
        factor(levels = c("WT", "HET", "HOM")),
      tg_parent = case_when(
        cross == "WT" ~ "Wildtype cross",
        str_detect(cross, "Female") ~ "Transgenic mother",
        str_detect(cross, "Male") ~ "Transgenic father"
      ) |>
        factor(
          levels = c("Transgenic mother", "Wildtype cross", "Transgenic father")
        )
    )
}

# Genotype colours; same values as `treatment_colours` in
# blood_feeding_functions.R so Figures 1 and 2 match
genotype_colours <- c("WT" = "#4D4D4D", "HET" = "#E69F00", "HOM" = "#56B4E9")

parent_shapes <- c(
  "Transgenic mother" = 16,
  "Wildtype cross" = 15,
  "Transgenic father" = 17
)

# Parsed facet labels; panels ordered so the two site-384 lines sit together
line_labels_parsed <- c(
  "1759" = "italic(cd)^g384",
  "QA383P" = "italic(cd)^{384*R}",
  "2360B5" = "italic(cd)^g225"
)
facet_line_order <- names(line_labels_parsed)

#' Plot a fertility outcome faceted by line (Figure 1)
#'
#' One panel per line, genotype on the x axis (coloured as in Figure 2),
#' cross direction shown by shape, and a dashed reference line at that line's
#' own wildtype estimate. Model estimates (emmeans, asymptotic 95% CI) are
#' overlaid on jittered raw observations.
#'
#' @param emm_data A tibble from [get_emm()] with columns `line`, `cross`,
#'   `response`, `asymp.LCL`, `asymp.UCL`.
#' @param raw_data A data frame of raw observations with `line` and `cross`.
#' @param raw_y <tidy-select> Column in `raw_data` for the raw points.
#' @param y_label Character. Y-axis label.
#' @param y_scale A ggplot2 y scale, or `NULL`. Default `NULL`.
#' @return A `ggplot` object.
plot_fertility_by_line <- function(
  emm_data,
  raw_data,
  raw_y,
  y_label,
  y_scale = NULL
) {
  relevel_line <- \(df) {
    mutate(df, line = factor(as.character(line), levels = facet_line_order))
  }
  emm_data <- emm_data |> split_cross() |> relevel_line()
  raw_data <- raw_data |> split_cross() |> relevel_line()
  wt_ref <- filter(emm_data, genotype == "WT")
  dodge <- position_dodge(width = 0.6)

  p <- ggplot(
    emm_data,
    aes(
      x = genotype,
      y = response,
      colour = genotype,
      shape = tg_parent,
      group = tg_parent
    )
  ) +
    geom_hline(
      data = wt_ref,
      aes(yintercept = response),
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.4
    ) +
    geom_point(
      data = raw_data,
      aes(y = {{ raw_y }}),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.6,
        seed = 1
      ),
      alpha = 0.25,
      size = 1.2
    ) +
    geom_errorbar(
      aes(ymin = asymp.LCL, ymax = asymp.UCL),
      position = dodge,
      width = 0.25,
      linewidth = 0.8
    ) +
    geom_point(position = dodge, size = 3) +
    facet_wrap(
      ~line,
      labeller = as_labeller(line_labels_parsed, label_parsed),
      drop = FALSE
    ) +
    scale_colour_manual(values = genotype_colours, name = "Genotype") +
    scale_shape_manual(values = parent_shapes, name = "Cross") +
    guides(
      colour = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    labs(x = "Genotype", y = y_label) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  if (!is.null(y_scale)) {
    p <- p + y_scale
  }
  p
}
