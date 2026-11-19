# Transgenic insertions in cardinal cause fitness costs in *Anopheles stephensi*

**Analytical code for the *cardinal* fitness manuscript**

> Authors: Mireia Larrosa-Godall, Lewis Shackleford, Philip T. Leftwich,
> Estela Gonzalez, Joshua X. D. Ang, Matthew Edwards, Katherine Nevard,
> James C.Y. Luk, Morgan Mckee, Rob Noad, Michelle A. E. Anderson, Luke Alphey
>
> Correspondence: michelle.anderson@york.ac.uk, luke.alphey@york.ac.uk

---

## Licence

This code is released under the [MIT License](LICENSE). You are free to reuse
and adapt it with attribution.

---

## Associated Manuscript

Larrosa-Godall, M. *et al.* (in prep). Transgenic insertions in cardinal cause fitness costs in *Anopheles stephensi*.

---

## Software

All analyses were run under **R v4.5.2**, with package versions pinned in
`renv.lock`. Rendering the report also needs [Quarto](https://quarto.org/).

Key packages (see `scripts/packages.R` for the loaded set; `gt` is loaded
by `results_summary.qmd` and `scripts/fitness_summary_table.R`):

| Package   | Version (`renv.lock`) | Role                                            |
|-----------|-----------------------|-------------------------------------------------|
| tidyverse | 2.0.0                 | Data manipulation and visualisation             |
| here      | 1.0.2                 | Portable file paths                             |
| readxl    | 1.5.0                 | Excel file import                               |
| janitor   | 2.2.1                 | Data cleaning utilities                         |
| glmmTMB   | 1.1.14                | Zero-inflated and mixed models                  |
| flexsurv  | 2.3.2                 | Parametric (Weibull) survival models            |
| survival  | 3.8-3                 | Kaplan–Meier estimates and log-rank tests       |
| emmeans   | 2.0.3                 | Estimated marginal means and contrasts          |
| ggplot2   | 4.0.3                 | Data visualisation                              |
| patchwork | 1.3.2                 | Plot composition                                |
| assertr   | 3.0.1                 | Data validation                                 |
| gt        | 1.3.0                 | Report and manuscript tables                    |
| sjPlot    | 2.9.0                 | Model coefficient tables                        |

> **Before running:** restore the package library with `renv::restore()`
> (`.Rprofile` activates `renv` automatically when the project is opened).

---

## Data

Data are stored in `data/`. The scripts read these files; the only file a
script writes into `data/` is `ddCt_results.rds` (see below).

| File | Contents | Used by |
|------|----------|---------|
| `FertilityData.xlsx` | Fecundity and hatching-rate data | `fertility_data.R` |
| `SurvivalAfterBFData.xlsx` | Survival after blood feeding | `blood_feeding.R` |
| `TransHetSurvivalData.xlsx` | Trans-heterozygote survival | `blood_feeding.R` |
| `SurvivalData.xlsx` | Larva-to-adult survival (eclosion) | `hom_survival.R` |
| `smurf.xlsx` | Smurf assay counts | `smurf.R` |
| `SmurfAssayData.xlsx` | Survival after 'blue-feeding' and by Smurf status | `smurf.R` |
| `XASupplementationData.xlsx` | Xanthurenic acid supplementation survival | `supplementation_data.R` |
| `Cardinal_Rps7_Results_20250611.csv` | qPCR Cq values, rps7 reference | `delta_delta_Ct_script.R` |
| `Cardinal_GADPH_Results_20250613.csv` | qPCR Cq values, GAPDH reference | `delta_delta_Ct_script.R` |
| `qPCR.xlsx` | qPCR workbook (not read by any script; see the CSV files) | none |
| `ddCt_results.rds` | **Output** of `delta_delta_Ct_script.R`: qPCR plots and contrasts, read by the report and `save_figures.R` | `results_summary.qmd`, `save_figures.R` |



---

## Repository Structure

```
cardinal_fitness/
├── data/                              # Input data (and ddCt_results.rds output)
|
├── scripts/
│   ├── packages.R                     # All package loading (sourced by every script)
│   ├── fertility_data.R               # Fecundity and hatching rate (Figure 1)
│   ├── hom_survival.R                 # Larva-to-adult survival (eclosion)
│   ├── blood_feeding.R                # Survival after blood feeding (Figure 2) and trans-hets (Figure 5)
│   ├── smurf.R                        # Smurf assay and survival (Figure 3, SFig. 3)
│   ├── supplementation_data.R         # XA supplementation (Figure 4)
│   ├── delta_delta_Ct_script.R        # qPCR delta-delta-Ct (Figure 6, SFig. 5)
│   ├── fitness_summary_table.R        # Table 1 and the supplementary fitness table
│   ├── save_figures.R                 # Exports all figures to outputs/figures/
│   └── functions/
│       ├── fertility_data_functions.R
│       ├── blood_feeding_functions.R
│       ├── smurf_functions.R
│       ├── supplementation_data_functions.R
│       ├── delta_delta_Ct_script_functions.R
│       └── survival_helpers.R         # Shared helpers for the survival analyses
├── results_summary.qmd                # Quarto report of all model tables
├── results_summary.html               # Rendered report
├── renv.lock                          # Package versions
├── air.toml                           # Code formatter settings
├── project.Rproj
├── LICENSE
└── README.md
```



---

## Run Order

Each analysis script sources `scripts/packages.R` and its own functions
file(s). The analysis scripts are independent of one another, with these
exceptions:

- `delta_delta_Ct_script.R` must be run before `save_figures.R` and
  `results_summary.qmd`, because it creates `data/ddCt_results.rds`.
- `fitness_summary_table.R` sources `fertility_data.R`, `hom_survival.R`,
  `blood_feeding.R` and `smurf.R`.
- `save_figures.R` sources `fertility_data.R`, `blood_feeding.R`, `smurf.R`
  and `supplementation_data.R`.
- `results_summary.qmd` sources `fertility_data.R`, `blood_feeding.R`,
  `smurf.R`, `supplementation_data.R`, `hom_survival.R` and
  `fitness_summary_table.R`.

Recommended order:

1. `scripts/delta_delta_Ct_script.R`: **Figure 6** and **SFig. 5** (qPCR, both reference genes); saves `data/ddCt_results.rds`
2. `scripts/fitness_summary_table.R`: **Table 1** and the supplementary fitness table (`outputs/tables/`)
3. `scripts/save_figures.R`: **Figures 1–6** and **SFig. 5** (`outputs/figures/`)
4. Render the report (below)

Individual figures come from: `fertility_data.R` (Figure 1), `blood_feeding.R`
(Figures 2 and 5), `smurf.R` (Figure 3 and SFig. 3), `supplementation_data.R`
(Figure 4). Figure 3 photographs and the feeding-rate table (panels A and B)
are assembled outside R.


**Computationally intensive steps:**

- `delta_delta_Ct_script.R` runs 5,000 permutations per reference-gene
  analysis (Freedman–Lane permutation within sex; seed set in the script).
  Expect a few minutes per run on a modern laptop.
- The distribution-comparison block in `supplementation_data.R` fits several
  parametric models with `flexsurvreg`, which takes a few minutes.

---

## Inputs and Outputs

| Script | Input(s) | Main output(s) |
|--------|----------|----------------|
| `fertility_data.R` | `FertilityData.xlsx` | `p_fecundity`, `p_fertility`, `fig1` (ggplot objects); `egg_model`, `larvae_model2`; `egg_vs_wt`, `larvae_vs_wt` (contrasts against each line's WT) |
| `hom_survival.R` | `SurvivalData.xlsx` | `model`; `eclosion_observed`; `eclosion_vs_wt` |
| `blood_feeding.R` | `SurvivalAfterBFData.xlsx`, `TransHetSurvivalData.xlsx` | `survival_after_blood`, `transhet_survival` (ggplot objects); `weibull_bf`, `weibull_th`; `bf_vs_wt`, `th_contrasts` |
| `smurf.R` | `smurf.xlsx`, `SmurfAssayData.xlsx` | `p_smurf`, `survival_after_smurf_plot`, `survival_by_smurf_plot` (ggplot objects); `smurf_model`; `model_smurf_surv`, `model_smurf_by` |
| `supplementation_data.R` | `XASupplementationData.xlsx` | `XA_suppl` (ggplot object); `xa_observed` and per-route contrast tables |
| `delta_delta_Ct_script.R` | `Cardinal_Rps7_Results_20250611.csv`, `Cardinal_GADPH_Results_20250613.csv` | `results` list (plots and contrasts), saved as `data/ddCt_results.rds` |
| `fitness_summary_table.R` | the objects above | `outputs/tables/table1_fitness_summary.{html,docx}`, `outputs/tables/tableS_fitness_vs_wt.{html,docx,csv}` |
| `save_figures.R` | the objects above and `data/ddCt_results.rds` | `outputs/figures/fig1`–`fig6` and `sfig5` (`.png` and `.pdf`) |
| `results_summary.qmd` | the scripts above | `results_summary.html` |
