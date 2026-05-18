# Changelog

All notable changes to this project are documented here.
This project does not follow Semantic Versioning; entries are dated.

---

## [Unreleased]

### Fixed
- Moved all function-definition files from `scripts/` to `scripts/functions/`
  so that function definitions and analysis scripts are cleanly separated
  (code review finding, Section 7).
- Replaced deprecated `map_df()` with `map() |> list_rbind()` in
  `supplementation_data.R` (code review finding, Section 5f).
- Replaced `group_by()` / `ungroup()` blocks with `.by` argument inside
  `summarise()` / `mutate()` in `smurf.R` and
  `delta_delta_Ct_script_functions.R` (code review finding, Section 5d).
- Replaced character-vector join syntax with `join_by()` in
  `delta_delta_Ct_script_functions.R` (code review finding, Section 5c).
- Replaced superseded `separate()` calls with `separate_wider_delim()` across
  multiple scripts (code review finding, Section 5h).
- Fixed duplicate `experiment` term in the `model_4` formula in
  `supplementation_data.R` (code review finding, Section 8).
- Fixed missing space before `+` on `transhet_theme+` line in
  `blood_feeding_functions.R` (code review finding, Section 5i).
- Added `here::here()` to bare `source()` call in `hom_survival.R`; replaced
  bare relative data paths with `here::here()` across multiple scripts
  (code review finding, Section 2).

### Added
- MIT `LICENSE` file (code review finding, Section 3 — blocking).
- Proper `README.md` with author contact, manuscript title, R version, key
  packages, run order, computationally intensive step warnings, and
  input/output descriptions (code review finding, Section 3).
- `CHANGELOG.md` (this file) (code review finding, Section 1).
- `outputs/figures/` and `outputs/tables/` directories (code review finding,
  Section 1).
- `assertr` biological plausibility validation blocks in data-loading scripts
  (code review finding, Section 9 — blocking).
- Added `library(assertr)` to `scripts/packages.R`.
- Added roxygen2 documentation block to `get_emm_contrasts()` in
  `fertility_data_functions.R` (code review finding, Section 7).

### Changed
- All analysis scripts updated to source function files from
  `scripts/functions/` (accompanies function-file move above).

---

## Notes

- `renv.lock` is absent from the repository. Before manuscript submission,
  run `renv::init()` and `renv::snapshot()` to generate a lockfile and
  commit it. Collaborators should run `renv::restore()` before executing
  the pipeline.
