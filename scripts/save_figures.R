# save_figures.R
# Export the manuscript figures produced by the analysis scripts to
# outputs/figures/ as PDF (vector, cairo for Unicode legend text) and PNG.
# Figure 3 photographs and the feeding-rate table (panels A/B) are assembled
# outside R; only its R-generated panels are saved here.

source(here::here("scripts", "packages.R"))

fig_dir <- here::here("outputs", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(plot, name, width, height) {
  ggsave(
    file.path(fig_dir, paste0(name, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    file.path(fig_dir, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# Figure 1: fecundity (A) and fertility (B), faceted by line
source(here::here("scripts", "fertility_data.R"))
save_fig(fig1, "fig1_fecundity_fertility", width = 10, height = 8)

# Figure 2: survival after blood feeding
# Figure 5: trans-heterozygote survival
source(here::here("scripts", "blood_feeding.R"))
save_fig(survival_after_blood, "fig2_survival_after_blood_feeding", width = 9, height = 6)
save_fig(transhet_survival, "fig5_transhet_survival", width = 10, height = 4.5)

# Figure 3: Smurf proportions and survival after 'blue-feeding'
source(here::here("scripts", "smurf.R"))
save_fig(p_smurf, "fig3_smurf_proportion", width = 5, height = 4)
save_fig(survival_after_smurf_plot, "fig3_survival_by_genotype", width = 9, height = 4)
save_fig(survival_by_smurf_plot, "fig3_survival_by_smurf_status", width = 6, height = 4)

# Figure 4: XA supplementation
source(here::here("scripts", "supplementation_data.R"))
save_fig(XA_suppl, "fig4_xa_supplementation", width = 9, height = 6)

# Figure 6 (cd / rps7) and Supplementary Figure 5 (cd / GAPDH): qPCR, from the
# results saved by scripts/delta_delta_Ct_script.R (corrected Freedman-Lane
# permutation p-values)
qpcr <- readRDS(here::here("data", "ddCt_results.rds"))
save_fig(qpcr$Cardinal_Rps7$plot, "fig6_qpcr_rps7", width = 8, height = 5)
save_fig(qpcr$Cardinal_GADPH$plot, "sfig5_qpcr_gapdh", width = 8, height = 5)
