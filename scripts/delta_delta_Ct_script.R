source(here::here("scripts", "packages.R"))

# ── Shared lookups ─────────────────────────────────────────────────────────────

sample_levels <- c("2360B5", "1759", "D251", "QA383P", "SDA-500")

sample_name <- c(
  "2360B5" = expression(italic(cd)^g225),
  "1759" = expression(italic(cd)^g384),
  "D251" = expression(italic(cd)^225),
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

# Full delta-delta-Ct pipeline: raw CSV → tidy combined_data with
# Relative_Expression, Sample (factored), Sex columns
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

# Standard ggplot for relative expression
plot_ddct <- function(data, contrast, plot_title = NULL) {
  ggplot(data, aes(x = Sample, y = delta_delta_Ct, fill = Sample)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_grid(. ~ Sex, scales = "free_x", space = "free_x") +
    stat_pvalue_manual(contrast) +
    scale_x_discrete(labels = sample_name) +
    scale_fill_manual(values = id_colors) +
    labs(
      x = NULL,
      y = expression(Delta * Delta * italic(C)[T]),
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

# GLM + emmeans treatment-vs-control contrasts
contrast_ddct <- function(data, control = "SDA-500") {
  model <- glm(delta_delta_Ct ~ Sample * Sex, data = data, family = gaussian)
  emm <- emmeans(model, ~ Sample | Sex)
  contrast(emm, method = "trt.vs.ctrl", ref = control, adjust = "none")
}

# Wrangle emmeans contrast output into stat_pvalue_manual format
tidy_contrast <- function(contrasts, y_positions) {
  summary(contrasts) |>
    as_tibble() |>
    mutate(
      p.value = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        .default = "n.s."
      )
    ) |>
    separate(contrast, into = c("group1", "group2"), sep = " - ") |>
    mutate(
      group2 = str_remove_all(group2, "[()]"),
      y.position = y_positions[group1]
    )
}

# Residual diagnostics for the GLM underlying contrast_ddct
check_ddct <- function(data) {
  model <- glm(delta_delta_Ct ~ Sample * Sex, data = data, family = gaussian)
  par(mfrow = c(2, 2))
  plot(model)
  par(mfrow = c(1, 1))
  model
}

# ── Analysis ──────────────────────────────────────────────────────────────────

# Input specifications for each reference gene
analyses <- list(
  Cardinal_Rps7 = list(
    file = "data/Cardinal_Rps7_Results_20250611.csv",
    ref_gene = "rps7",
    title = "Cardinal / rps7",
    y_pos = c("2360B5" = 6, "1759" = 7, "D251" = 8, "QA383P" = 9)
  ),
  Cardinal_GADPH = list(
    file = "data/Cardinal_GADPH_Results_20250613.csv",
    ref_gene = "GAPDH",
    title = "Cardinal / GAPDH",
    y_pos = c("2360B5" = 8, "1759" = 9, "D251" = 10, "QA383P" = 11)
  )
)

# Run pipeline for each analysis
results <- map(analyses, function(a) {
  data <- calc_ddct(a$file, target_gene = "Cardinal", ref_gene = a$ref_gene)
  contrasts <- contrast_ddct(data)
  contrast <- tidy_contrast(contrasts, a$y_pos)
  plot <- plot_ddct(data, contrast, plot_title = a$title)
  list(data = data, contrast = contrast, plot = plot)
})

# Access results
results$Cardinal_Rps7$plot
results$Cardinal_GADPH$plot

# ── Residual checks ───────────────────────────────────────────────────────────

check_ddct(results$Cardinal_Rps7$data) # Cardinal / rps7
check_ddct(results$Cardinal_GADPH$data) # Cardinal / GAPDH
