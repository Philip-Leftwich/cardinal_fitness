source(here::here("scripts", "packages.R"))
source(here::here("scripts", "delta_delta_Ct_script_functions.R"))

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
