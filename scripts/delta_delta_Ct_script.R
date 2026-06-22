source(here::here("scripts", "packages.R"))
source(here::here("scripts", "functions", "delta_delta_Ct_script_functions.R"))

# ── Analysis ──────────────────────────────────────────────────────────────────
set.seed(220626) # For reproducibility of permutation tests
# Input specifications for each reference gene
analyses <- list(
  Cardinal_Rps7 = list(
    file = here::here("data", "Cardinal_Rps7_Results_20250611.csv"),
    ref_gene = "rps7",
    title = "Cardinal / rps7",
    y_pos = c("2360B5" = 9, "1759" = 8, "D251" = 7, "QA383P" = 6)
  ),
  Cardinal_GADPH = list(
    file = here::here("data", "Cardinal_GADPH_Results_20250613.csv"),
    ref_gene = "GAPDH",
    title = "Cardinal / GAPDH",
    y_pos = c("2360B5" = 9, "1759" = 8, "D251" = 7, "QA383P" = 6)
  )
)

# Run pipeline for each analysis
results <- map(analyses, function(a) {
  data <- calc_ddct(a$file, target_gene = "Cardinal", ref_gene = a$ref_gene)
  contrasts <- contrast_ddct(data)
  contrast <- tidy_contrast(contrasts$permuted, a$y_pos)
  plot <- plot_ddct(data, contrast, plot_title = a$title)
  list(data = data, contrasts = contrasts, contrast = contrast, plot = plot)
})

# Access results
results$Cardinal_Rps7$plot
results$Cardinal_GADPH$plot

saveRDS(results, here::here("data", "ddCt_results.rds"))

# ── Residual checks ───────────────────────────────────────────────────────────

check_ddct(results$Cardinal_Rps7$data) # Cardinal / rps7
check_ddct(results$Cardinal_GADPH$data) # Cardinal / GAPDH
