source(here::here("scripts", "packages.R"))

# ── Shared lookups ─────────────────────────────────────────────────────────────

sample_levels <- c("SDA-500", "2360B5", "1759", "D251", "QA383P")

sample_name <- c(
  "2360B5" = expression(italic(cd)^g225),
  "1759" = expression(italic(cd)^g384),
  "D251" = expression(
    italic(cd)^{
      "225R"
    }
  ),
  "QA383P" = expression(
    italic(cd)^{
      "384R"
    }
  ),
  "SDA-500" = expression(WT)
)

id_colors <- c(
  "2360B5" = "#8BABD3",
  "1759" = "#BB8BD3",
  "D251" = "#FFA040",
  "QA383P" = "#FFE699",
  "SDA-500" = "lightgray"
)

# ── Helper functions ───────────────────────────────────────────────────────────

#' Run the full delta-delta-Ct pipeline
#'
#' Reads a raw qPCR CSV, averages technical replicates, computes delta-Ct
#' (target minus reference), then delta-delta-Ct relative to a control sample,
#' and returns 2^(-ddCt) as relative expression.
#'
#' @param file Path to a CSV file with columns `Sample`, `Target`, `Replicate`,
#'   and `Cq`.
#' @param target_gene Character. Name of the target gene matching values in the
#'   `Target` column.
#' @param ref_gene Character. Name of the reference/housekeeping gene.
#' @param control_sample Character. The sample label used as the ddCt baseline.
#'   Default `"SDA-500 fem"`.
#' @return A tibble with columns including `Sample`, `Sex`, `delta_delta_Ct`,
#'   and `Relative_Expression`, with `Sample` factored to `sample_levels`.
calc_ddct <- function(
  file,
  target_gene,
  ref_gene,
  control_sample = "SDA-500 fem"
) {
  raw <- read.csv(file) |>
    filter(Sample != "NTC") |>
    mutate(Cq = as.numeric(Cq)) |>
    group_by(Sample, Target, Replicate) |>
    summarise(Cq.mean = mean(Cq), .groups = "drop")

  target <- raw |> filter(Target == target_gene)
  ref <- raw |> filter(Target == ref_gene) |> rename(Cq.ref = Cq.mean)

  left_join(target, ref, by = c("Sample", "Replicate")) |>
    mutate(delta_Ct = Cq.mean - Cq.ref) |>
    group_by(Sample, Replicate) |>
    mutate(mean_deltaCt = mean(delta_Ct)) |>
    ungroup() |>
    mutate(
      mean_control = mean(mean_deltaCt[Sample == control_sample]),
      delta_delta_Ct = delta_Ct - mean_control,
      Relative_Expression = 2^(-delta_delta_Ct)
    ) |>
    separate(Sample, into = c("Sample", "Sex"), sep = " ", remove = TRUE) |>
    mutate(
      Sex = recode(Sex, "fem" = "Female", "female" = "Female", "male" = "Male"),
      Sample = factor(Sample, levels = sample_levels)
    )
}

#' Plot delta-delta-Ct relative expression
#'
#' Produces a faceted boxplot (by `Sex`) of delta-delta-Ct values with
#' overlaid jitter and p-value annotations.
#'
#' @param data A tibble from [calc_ddct()] containing `Sample`, `Sex`, and
#'   `delta_delta_Ct` columns.
#' @param contrast A tibble from [tidy_contrast()] formatted for
#'   [ggpubr::stat_pvalue_manual()].
#' @param plot_title Character or `NULL`. Optional plot title. Default `NULL`.
#' @return A `ggplot` object.
plot_ddct <- function(data, contrast, plot_title = NULL) {
  ggplot(data, aes(x = Sample, y = Relative_Expression, fill = Sample)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_grid(. ~ Sex, scales = "free_y") +
    stat_pvalue_manual(contrast) +
    scale_y_continuous(
      trans = "log2",
      name = expression(
        log[2] ~ "Relative Expression (" *
          2^{
            -Delta * Delta * "Ct"
          } *
          ")"
      ),
      breaks = scales::breaks_log(n = 6, base = 2),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_x_discrete(labels = sample_name, limits = rev(levels(data$Sample))) +
    scale_fill_manual(values = id_colors) +
    labs(
      x = NULL,
      title = plot_title
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 12),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "none"
    )
}

#' Compute treatment-vs-control contrasts from a ddCt GLM
#'
#' Fits a Gaussian GLM of `delta_delta_Ct ~ Sample * Sex` and returns
#' pairwise treatment-vs-control contrasts via [emmeans::contrast()].
#'
#' @param data A tibble from [calc_ddct()].
#' @param control Character. The reference level in `Sample` to contrast
#'   against. Default `"SDA-500"`.
#' @return An `emmGrid` contrast object from [emmeans::contrast()].
contrast_ddct <- function(data, control = "SDA-500", nperms = 5000) {
  # Full model
  model_full <- lm(delta_Ct ~ Sample * Sex, data = data)
  emm <- emmeans(model_full, ~ Sample | Sex)
  con_obs <- contrast(emm, method = "trt.vs.ctrl", adjust = "none")
  obs_stats <- summary(con_obs)$t.ratio

  # Reduced model — removes Sample effect, keeps Sex
  model_reduced <- lm(delta_Ct ~ Sex, data = data)

  permute <- function(data, model_full, model_reduced) {
    weights <- sample(c(-1, 1), nrow(data), replace = TRUE)
    data$delta_Ct_perm <- fitted(model_reduced) +
      residuals(model_full) * weights
    fit_perm <- lm(delta_Ct_perm ~ Sample * Sex, data = data)
    emm_perm <- emmeans(fit_perm, ~ Sample | Sex)
    con_perm <- contrast(emm_perm, method = "trt.vs.ctrl", adjust = "none")
    summary(con_perm)$t.ratio
  }

  set.seed(42)
  perm_dist <- replicate(nperms, permute(data, model_full, model_reduced))
  p_raw <- rowMeans(abs(perm_dist) >= abs(obs_stats))

  con_parametric <- summary(con_obs) |>
    as_tibble()

  con_permuted <- con_parametric |>
    mutate(p.perm = p_raw)

  diagnostics <- tibble(
    contrast = con_parametric$contrast,
    sex = con_parametric$Sex,
    obs_t = obs_stats,
    perm_mean = rowMeans(perm_dist),
    perm_sd = apply(perm_dist, 1, sd)
  )

  results <- list(
    parametric = con_parametric,
    permuted = con_permuted,
    diagnostics = diagnostics
  )

  return(results)
}

#' Tidy emmeans contrasts for use with stat_pvalue_manual
#'
#' Converts an emmeans contrast summary into a tibble formatted for
#' [ggpubr::stat_pvalue_manual()], with significance stars and y-axis
#' positions.
#'
#' @param contrasts An `emmGrid` contrast object from [contrast_ddct()].
#' @param y_positions A named numeric vector mapping `group1` sample names to
#'   y-axis positions for the bracket annotations.
#' @return A tibble with columns `group1`, `group2`, `p.value` (as
#'   significance stars), and `y.position`.
tidy_contrast <- function(contrasts, y_positions) {
  report_p <- function(p, digits = 3) {
    reported <- if_else(p < 0.001, "p < 0.001", paste("p =", round(p, digits)))

    return(reported)
  }

  contrasts |>
    mutate(
      p.value = report_p(p.perm)
    ) |>
    separate(contrast, into = c("group1", "group2"), sep = " - ") |>
    mutate(
      group2 = str_remove_all(group2, "[()]"),
      y.position = y_positions[group1]
    )
}

#' Plot residual diagnostics for the ddCt GLM
#'
#' Fits the same Gaussian GLM as [contrast_ddct()] and displays a 2×2
#' base-R diagnostic plot (residuals vs fitted, Q-Q, scale-location,
#' leverage).
#'
#' @param data A tibble from [calc_ddct()].
#' @return The fitted `glm` object, invisibly.
check_ddct <- function(data) {
  model <- glm(delta_Ct ~ Sample * Sex, data = data, family = gaussian)
  par(mfrow = c(2, 2))
  plot(model)
  par(mfrow = c(1, 1))
  model
}
