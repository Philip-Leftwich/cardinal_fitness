# cardinal_fitness

**Analytical code for the *cardinal* fitness manuscript**

> Authors: Mireia Larrosa-Godall, Lewis Shackleford, Philip T. Leftwich,
> Estela Gonzalez, Joshua X. D. Ang, Matthew Edwards, Katherine Nevard,
> James C.Y. Luk, Morgan Mckee, Michelle A. E. Anderson, Luke Alphey
>
> Correspondence: m.anderson@pirbright.ac.uk; l.alphey@pirbright.ac.uk



---

## Licence

This code is released under the [MIT License](LICENSE). You are free to reuse
and adapt it with attribution.

---

## Associated Manuscript

Larrosa-Godall, M. *et al.* (in prep). *Fitness effects of cardinal
mutations in* Anopheles stephensi. (Provisional title —
update before submission.)


---

## R Version and Key Packages

All analyses were run under **R v4.5.2**.

Key packages (see `scripts/packages.R` for the full list):

| Package    | Version | Role                                   |
|------------|---------|----------------------------------------|
| tidyverse  | 2.0.0   | Data manipulation and visualisation    |
| here       | 1.0.1   | Portable file paths                    |
| readxl     | 1.4.3   | Excel file import                      |
| janitor    | 2.2.0   | Data cleaning utilities                |
| glmmTMB    | 1.1.9   | Zero-inflated and mixed models         |
| flexsurv   | 2.3     | Parametric survival models             |
| survival   | 3.5-8   | Kaplan–Meier and Cox models            |
| emmeans    | 1.10.1  | Estimated marginal means and contrasts |
| ggplot2    | 3.5.1   | Data visualisation                     |
| patchwork  | 1.2.0   | Plot composition                       |
| assertr    | 2.9.0   | Data validation                        |

Full session information is printed at the end of `report.qmd`.

> **Before running:** restore the package library with `renv::restore()`.
> A `renv.lock` file should be committed with the submission snapshot.
> If it is absent, run `renv::init()` followed by `renv::snapshot()` to
> generate it.

---

## Data

Raw data are stored in `data/` (Excel and CSV files). All data files are
read by the scripts; no script writes to `data/raw/`.

| File | Contents |
|------|----------|
| `FertilityData.xlsx` | Fecundity and hatching-rate data |
| `SurvivalAfterBFData.xlsx` | Survival after blood feeding |
| `TransHetSurvivalData.xlsx` | Trans-heterozygote survival |
| `SurvivalData.xlsx` | Homozygous viability (eclosion) |
| `SmurfAssayData.xlsx` | Smurf midgut-integrity assay |
| `XASupplementationData.xlsx` | Xanthurenic acid supplementation survival |
| `smurf.xlsx` | Smurf screening raw data |
| `qPCR.xlsx` | qPCR raw data (unused by pipeline; see CSV files) |
| `Cardinal_Rps7_Results_20250611.csv` | qPCR Cq values — rps7 reference |
| `Cardinal_GADPH_Results_20250613.csv` | qPCR Cq values — GAPDH reference |

Data are not archived separately at this time. A persistent identifier
(DOI) should be added here before manuscript submission.

---

## Repository Structure

```
cardinal_fitness/
├── data/                          # Raw data files (read-only in pipeline)
├── outputs/
│   ├── figures/                   # Saved figures (.png / .pdf)
│   └── tables/                    # Saved model-summary tables
├── scripts/
│   ├── packages.R                 # All package loading (source first)
│   ├── blood_feeding.R            # Survival after blood feeding
│   ├── delta_delta_Ct_script.R    # qPCR delta-delta-Ct analysis
│   ├── fertility_data.R           # Fecundity and hatching-rate analysis
│   ├── hom_survival.R             # Homozygous viability (eclosion) analysis
│   ├── smurf.R                    # Smurf assay and smurf-type survival
│   ├── supplementation_data.R     # XA supplementation survival
│   └── functions/
│       ├── blood_feeding_functions.R
│       ├── delta_delta_Ct_script_functions.R
│       ├── fertility_data_functions.R
│       ├── smurf_functions.R
│       └── supplementation_data_functions.R
├── report.qmd                     # Quarto analysis report (renders HTML)
├── CHANGELOG.md
├── LICENSE
└── README.md
```

---

## Run Order

Scripts are independent (no script depends on the output of another,
except `delta_delta_Ct_script.R` which saves `data/ddCt_results.rds`).
Each analysis script sources `scripts/packages.R` and its own functions
file automatically.

The recommended run order (matching the figures in the manuscript):

1. `scripts/fertility_data.R` — **Figure 1** (fecundity and hatching rate)
2. `scripts/blood_feeding.R` — **Figure 2** (survival after blood feeding) and **Figure 5** (trans-het survival)
3. `scripts/smurf.R` — **Figure 3** (smurf assay and smurf-type survival)
4. `scripts/supplementation_data.R` — **Figure 4** (XA supplementation)
5. `scripts/delta_delta_Ct_script.R` — **Figure 5** (qPCR, both reference genes)
6. `scripts/hom_survival.R` — **Supplementary** (homozygous eclosion)

To render the full results report: open `report.qmd` in RStudio and click
**Render**, or run `quarto render report.qmd` from the project root.

**Computationally intensive steps:**

- `delta_delta_Ct_script.R` runs 5 000 permutations per reference-gene
  analysis. Expect ~2–5 minutes per run on a modern laptop.
- The distribution-comparison block in `supplementation_data.R` fits five
  parametric models via `flexsurvreg`; allow ~5 minutes.

---

## Inputs and Outputs

| Script | Input(s) | Output(s) |
|--------|----------|-----------|
| `fertility_data.R` | `FertilityData.xlsx` | `p_fecundity`, `p_fertility` ggplot objects |
| `blood_feeding.R` | `SurvivalAfterBFData.xlsx`, `TransHetSurvivalData.xlsx` | `survival_after_blood`, `transhet_survival` ggplot objects |
| `smurf.R` | `smurf.xlsx`, `SmurfAssayData.xlsx` | `p_smurf`, `survival_after_smurf_plot`, `survival_by_smurf_plot` ggplot objects |
| `supplementation_data.R` | `XASupplementationData.xlsx` | `XA_suppl` ggplot object |
| `delta_delta_Ct_script.R` | `Cardinal_Rps7_Results_20250611.csv`, `Cardinal_GADPH_Results_20250613.csv` | `results` list (plots + contrasts); `data/ddCt_results.rds` |
| `hom_survival.R` | `SurvivalData.xlsx` | Model summary (printed to console) |
