# survival_helpers.R
# Shared helpers for the survival analyses and figures (Figs 2, 3, 4 and 5).
# Kept in a separate file so it can be sourced by any survival script without
# overwriting the script-specific lookups (e.g. `line_labels`, `tidy_km`)
# defined in the other function files.

# ── Estimate-type legend (Weibull model vs Kaplan-Meier) ───────────────────────

# Legend keys; map these as constant aesthetics, e.g.
#   geom_line(aes(linetype = lt_model)) and geom_step(aes(linetype = lt_km))
lt_model <- "Weibull model"
lt_km <- "Kaplan\u2013Meier"

estimate_linetypes <- setNames(c("solid", "dashed"), c(lt_model, lt_km))

#' Legend layers distinguishing Weibull model curves from Kaplan-Meier steps
#'
#' Replaces the old "Solid: Weibull model | Dashed: Kaplan-Meier" plot caption
#' with a linetype legend placed directly after the colour (genotype or
#' treatment) legend. The fill (CI ribbon) legend is suppressed because it
#' duplicates the colour legend and would split from it when a group has no
#' model ribbon (e.g. a group with no deaths shown by Kaplan-Meier only).
#'
#' @param colour_order Integer. Legend order for the colour legend. Default 1.
#' @param linetype_order Integer. Legend order for the estimate legend; set it
#'   higher than every colour legend in a patchwork with collected guides.
#'   Default 2.
#' @return A list of ggplot2 components to add to a plot.
estimate_linetype_layers <- function(colour_order = 1, linetype_order = 2) {
  list(
    scale_linetype_manual(
      values = estimate_linetypes,
      breaks = names(estimate_linetypes),
      name = "Estimate"
    ),
    guides(
      colour = guide_legend(
        order = colour_order,
        override.aes = list(linetype = "solid", linewidth = 0.8)
      ),
      fill = "none",
      linetype = guide_legend(
        order = linetype_order,
        override.aes = list(colour = "black", fill = NA, linewidth = 0.6)
      )
    )
  )
}

# ── Groups with no deaths (complete separation) ───────────────────────────────

#' Drop groups with no deaths before fitting a parametric survival model
#'
#' A group with no events has no finite maximum-likelihood time ratio
#' (complete separation). In `flexsurvreg` such a group also makes the
#' numerically estimated Hessian near-singular, which corrupts the standard
#' errors of *every* other coefficient (checked against
#' [survival::survreg()], which gives the correct standard errors). These
#' groups are therefore fitted out and reported descriptively as observed
#' counts; point estimates for the remaining groups are unchanged.
#'
#' @param data A data frame of individual survival records.
#' @param group <tidy-select> Grouping column (a factor or character).
#' @param event Character. Name of the event column (1 = died). Default
#'   `"event"`.
#' @return `data` without zero-event groups (unused factor levels dropped),
#'   with the removed group names in attribute `"zero_event_groups"`.
drop_zero_event_groups <- function(data, group, event = "event") {
  g <- rlang::as_name(rlang::enquo(group))
  deaths <- tapply(data[[event]], as.character(data[[g]]), sum)
  zero <- names(deaths)[deaths == 0]
  out <- data[!as.character(data[[g]]) %in% zero, ]
  if (is.factor(out[[g]])) {
    out[[g]] <- droplevels(out[[g]])
  }
  attr(out, "zero_event_groups") <- zero
  out
}

#' Observed sample size and deaths per group
#'
#' @param data A data frame of individual survival records.
#' @param group <tidy-select> Grouping column(s).
#' @param event Character. Name of the event column (1 = died, 0 = censored).
#'   Default `"event"`.
#' @return A tibble with columns for the group, `n` and `deaths`.
observed_deaths <- function(data, group, event = "event") {
  data |>
    summarise(n = n(), deaths = sum(.data[[event]]), .by = {{ group }})
}
