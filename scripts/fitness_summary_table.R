# fitness_summary_table.R
# Reviewer comment 5: summarise every fitness trait for each genotype as the
# direction and significance of its difference from the wildtype (WT) group of
# the same line, in one table (rows = genotype, columns = trait).
#
# Outputs (outputs/tables/):
#   table1_fitness_summary.{html,docx}   compact main-text table
#   tableS_fitness_vs_wt.{html,docx,csv} detailed supplementary table with
#                                        observed counts, ratio [95% CI], p
#
# All contrasts are against the same line's WT group, which differs by
# experiment: cd[g384] eclosion/fecundity = non-transgenic siblings;
# cd[384R] eclosion/fecundity = SDA-500; cd[g225] fecundity = siblings pooled
# with SDA-500 x SDA-500; survival after blood feeding = non-transgenic
# siblings of each line; Smurf assay = wildtype.

source(here::here("scripts", "packages.R"))
library(gt)

source(here::here("scripts", "fertility_data.R"))
source(here::here("scripts", "hom_survival.R"))
source(here::here("scripts", "blood_feeding.R"))
source(here::here("scripts", "smurf.R"))

out_dir <- here::here("outputs", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Row definitions ───────────────────────────────────────────────────────────
# `md` uses Pandoc superscripts (for .docx); `html` for the gt tables

rows <- tribble(
  ~row_key, ~group, ~md, ~html,
  "1759_HET", "Knock-in", "*cd*^g384^ HET", "<em>cd</em><sup>g384</sup> HET",
  "1759_HOM", "Knock-in", "*cd*^g384^ HOM", "<em>cd</em><sup>g384</sup> HOM",
  "2360B5_HET", "Knock-in", "*cd*^g225^ HET", "<em>cd</em><sup>g225</sup> HET",
  "2072_HET", "Knock-in", "*cd*^g384_del^ HET", "<em>cd</em><sup>g384_del</sup> HET",
  "2301_HET", "Knock-in", "*cd*^g338-384^ HET", "<em>cd</em><sup>g338-384</sup> HET",
  "QA383P_HET", "Knock-out", "*cd*^384R^ HET", "<em>cd</em><sup>384R</sup> HET",
  "QA383P_HOM", "Knock-out", "*cd*^384R^ HOM", "<em>cd</em><sup>384R</sup> HOM",
  "D251_HET", "Knock-out", "*cd*^225R^ HET", "<em>cd</em><sup>225R</sup> HET",
  "D251_HOM", "Knock-out", "*cd*^225R^ HOM", "<em>cd</em><sup>225R</sup> HOM",
  "TRANSHET", "Trans-heterozygote", "*cd*^g225^;*cd*^225R^", "<em>cd</em><sup>g225</sup>;<em>cd</em><sup>225R</sup>"
)

columns <- tribble(
  ~col_key, ~trait, ~sub, ~ratio_type,
  "eclosion_f", "Larva-to-adult survival", "Females", "OR",
  "eclosion_m", "Larva-to-adult survival", "Males", "OR",
  "eggs_mother", "Fecundity (eggs per female)", "Transgenic mother", "Rate ratio",
  "eggs_father", "Fecundity (eggs per female)", "Transgenic father", "Rate ratio",
  "hatch_mother", "Fertility (hatching rate)", "Transgenic mother", "OR",
  "hatch_father", "Fertility (hatching rate)", "Transgenic father", "OR",
  "survival", "Survival after blood meal", "Females", "Time ratio",
  "smurf", "Midgut barrier loss (Smurf+)", "Females", "OR"
)

fmt_n <- function(x, digits = 1) formatC(x, format = "f", digits = digits)
fmt_pct <- function(x) paste0(fmt_n(100 * x), "%")

# ── Larva-to-adult survival (eclosion) ────────────────────────────────────────

wt_eclosion <- eclosion_observed |>
  filter(Genotype == "WT") |>
  select(line, Sex, wt_adults = adults, wt_larvae = larvae)

res_eclosion <- eclosion_vs_wt |>
  mutate(Genotype = str_remove(contrast, " - WT")) |>
  left_join(eclosion_observed, by = c("line", "Genotype", "Sex")) |>
  left_join(wt_eclosion, by = c("line", "Sex")) |>
  transmute(
    row_key = paste(line, Genotype, sep = "_"),
    col_key = if_else(Sex == "female", "eclosion_f", "eclosion_m"),
    ratio,
    lcl = ratio.LCL,
    ucl = ratio.UCL,
    p = p.value,
    test = "Dunnett",
    observed = glue::glue(
      "{adults}/{larvae} adults ({fmt_pct(adults / larvae)}) vs WT ",
      "{wt_adults}/{wt_larvae} ({fmt_pct(wt_adults / wt_larvae)})"
    )
  )

# ── Fecundity and hatching ────────────────────────────────────────────────────

fert_obs <- fertility_data |>
  filter(!is.na(eggs)) |>
  summarise(
    n = n(),
    mean_eggs = mean(eggs),
    laying = mean(eggs > 0),
    .by = c(line, cross)
  )

hatch_obs <- fertility_data_combined |>
  filter(!is.na(larvae), eggs > 0) |>
  summarise(
    n_hatch = n(),
    larvae = sum(larvae),
    eggs = sum(eggs),
    .by = c(line, cross)
  )

fert_contrast_rows <- function(contrasts, kind) {
  # cross labels combine transgenic parent sex and genotype, e.g. "Female \nHET"
  contrasts |>
    mutate(
      cross = str_remove(contrast, " - WT$"),
      genotype = if_else(str_detect(cross, "HOM"), "HOM", "HET"),
      parent = if_else(str_detect(cross, "Female"), "mother", "father"),
      row_key = paste(line, genotype, sep = "_"),
      col_key = paste(kind, parent, sep = "_")
    )
}

wt_fert <- fert_obs |>
  filter(cross == "WT") |>
  select(line, wt_n = n, wt_mean = mean_eggs, wt_laying = laying)

res_eggs <- egg_vs_wt |>
  fert_contrast_rows("eggs") |>
  left_join(fert_obs, by = c("line", "cross")) |>
  left_join(wt_fert, by = "line") |>
  transmute(
    row_key,
    col_key,
    ratio,
    lcl = ratio.LCL,
    ucl = ratio.UCL,
    p = p.value,
    test = "Dunnett",
    observed = glue::glue(
      "mean {fmt_n(mean_eggs)} eggs, {fmt_pct(laying)} laid (n = {n}) vs WT ",
      "{fmt_n(wt_mean)}, {fmt_pct(wt_laying)} (n = {wt_n})"
    )
  )

wt_hatch <- hatch_obs |>
  filter(cross == "WT") |>
  select(line, wt_larvae = larvae, wt_eggs = eggs)

res_hatch <- larvae_vs_wt |>
  fert_contrast_rows("hatch") |>
  left_join(hatch_obs, by = c("line", "cross")) |>
  left_join(wt_hatch, by = "line") |>
  transmute(
    row_key,
    col_key,
    ratio,
    lcl = ratio.LCL,
    ucl = ratio.UCL,
    p = p.value,
    test = "Dunnett",
    observed = glue::glue(
      "{larvae}/{eggs} eggs hatched ({fmt_pct(larvae / eggs)}; ",
      "n = {n_hatch}) vs WT {fmt_pct(wt_larvae / wt_eggs)}"
    )
  )

# cd[g225] HET females died after blood feeding or did not feed (1 well, 0 eggs)
res_g225_mother <- tibble(
  row_key = "2360B5_HET",
  col_key = c("eggs_mother", "hatch_mother"),
  note = "not estimable: HET females died after blood feeding or did not feed"
)

# ── Survival after blood feeding ──────────────────────────────────────────────

# Log-rank test for comparisons where one group had no deaths (no finite
# time ratio)
logrank_p <- function(data, level, ref) {
  d <- filter(data, line_treatment %in% c(level, ref))
  fit <- survdiff(Surv(Hours, event) ~ line_treatment, data = d)
  pchisq(fit$chisq, df = 1, lower.tail = FALSE)
}

surv_observed <- function(deaths, n, deaths_ref, n_ref) {
  glue::glue("{deaths}/{n} died vs WT {deaths_ref}/{n_ref}")
}

res_bf <- bf_vs_wt |>
  mutate(
    row_key = str_replace(level, "_Het$", "_HET") |>
      str_replace("_Hom$", "_HOM"),
    separated = is.na(tr),
    p = if_else(
      separated,
      map2_dbl(level, ref, \(l, r) logrank_p(data_long, l, r)),
      p_adj
    )
  ) |>
  transmute(
    row_key,
    col_key = "survival",
    ratio = tr,
    lcl = tr_lcl,
    ucl = tr_ucl,
    p,
    test = if_else(separated, "log-rank", "Holm"),
    observed = surv_observed(deaths, n, deaths_ref, n_ref),
    note = case_when(
      separated ~ "no deaths in WT; time ratio not estimable",
      n < 5 ~ paste("only", n, "females"),
      .default = NA_character_
    )
  )

th_map <- c("D251:2360WT" = "D251_HET", "D251:2360Het" = "TRANSHET")

res_th <- th_contrasts |>
  filter(family == "vs_wt", level %in% names(th_map)) |>
  mutate(
    row_key = unname(th_map[level]),
    separated = is.na(tr),
    p = if_else(
      separated,
      map2_dbl(level, ref, \(l, r) logrank_p(data_transhet, l, r)),
      p_adj
    )
  ) |>
  transmute(
    row_key,
    col_key = "survival",
    ratio = tr,
    lcl = tr_lcl,
    ucl = tr_ucl,
    p,
    test = if_else(separated, "log-rank", "Holm"),
    observed = surv_observed(deaths, n, deaths_ref, n_ref),
    note = if_else(
      separated,
      "no deaths; time ratio not estimable (trans-het experiment; WT = cd[g225] siblings)",
      "trans-het experiment; WT = cd[g225] siblings"
    )
  )

# ── Smurf assay ───────────────────────────────────────────────────────────────

smurf_counts <- data_wide |>
  summarise(pos = sum(positive), neg = sum(negative), .by = genotype) |>
  mutate(n = pos + neg)
wt_s <- filter(smurf_counts, genotype == "WT")
ki_s <- filter(smurf_counts, genotype == "KI")
ko_s <- filter(smurf_counts, genotype == "KO")

ki_or <- broom::tidy(smurf_model, conf.int = TRUE, exponentiate = TRUE) |>
  filter(term == "genotypeKI")

ko_fisher <- fisher.test(matrix(c(ko_s$pos, ko_s$neg, wt_s$pos, wt_s$neg), 2, byrow = TRUE))

smurf_observed <- function(g) {
  glue::glue("{g$pos}/{g$n} Smurf+ ({fmt_pct(g$pos / g$n)}) vs WT {wt_s$pos}/{wt_s$n} ({fmt_pct(wt_s$pos / wt_s$n)})")
}

res_smurf <- tibble(
  row_key = c("2360B5_HET", "D251_HOM"),
  col_key = "smurf",
  ratio = c(ki_or$estimate, NA),
  lcl = c(ki_or$conf.low, NA),
  ucl = c(ki_or$conf.high, NA),
  p = c(ki_or$p.value, ko_fisher$p.value),
  test = c("Wald (GLM)", "Fisher's exact"),
  observed = c(smurf_observed(ki_s), smurf_observed(ko_s)),
  note = c(NA, "no Smurf+ females; odds ratio not estimable")
)

# ── Combine ───────────────────────────────────────────────────────────────────

fitness_results <- bind_rows(
  res_eclosion,
  res_eggs,
  res_hatch,
  res_g225_mother,
  res_bf,
  res_th,
  res_smurf
) |>
  mutate(observed = as.character(observed), note = as.character(note)) |>
  left_join(columns, by = "col_key") |>
  left_join(rows, by = "row_key")

stopifnot(!anyNA(fitness_results$trait), !anyNA(fitness_results$group)) # every cell mapped

# Cell for the compact table: direction relative to WT, ratio, significance.
# Direction is read from the ratio, or from the observed counts where the
# ratio is not estimable.
fitness_results <- fitness_results |>
  mutate(
    stars = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      .default = ""
    ),
    sep_dir = case_when(
      row_key == "2301_HET" & col_key == "survival" ~ -1, # 31/31 vs 0/47 died
      row_key == "D251_HET" & col_key == "survival" ~ 0,
      row_key == "D251_HOM" & col_key == "smurf" ~ 0,
      .default = NA
    ),
    direction = case_when(
      is.na(p) ~ NA_character_,
      p >= 0.05 ~ "≈",
      !is.na(ratio) & ratio < 1 ~ "↓",
      !is.na(ratio) & ratio > 1 ~ "↑",
      sep_dir < 0 ~ "↓",
      sep_dir > 0 ~ "↑"
    ),
    ratio_txt = if_else(
      is.na(ratio),
      "†",
      format(signif(ratio, 2), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
    ),
    flag = case_when(
      str_detect(coalesce(note, ""), "^only") ~ "§",
      .default = ""
    ),
    cell = case_when(
      is.na(p) ~ "n/e‡",
      .default = paste0(direction, " ", ratio_txt, stars, flag)
    )
  )

# ── Table 1 (compact) ─────────────────────────────────────────────────────────

wide <- fitness_results |>
  select(row_key, col_key, cell) |>
  pivot_wider(names_from = col_key, values_from = cell) |>
  right_join(rows, by = "row_key") |>
  mutate(row_key = factor(row_key, levels = rows$row_key)) |>
  arrange(row_key)

for (cc in columns$col_key) {
  if (!cc %in% names(wide)) wide[[cc]] <- NA_character_
}
wide <- wide |>
  mutate(across(all_of(columns$col_key), \(x) coalesce(x, "—")))

table1_footnote <- paste(
  "Each cell compares the genotype with the wildtype (WT) group of the same",
  "line and experiment: ↓/↑, significantly lower/higher than WT;",
  "≈, no significant difference (p ≥ 0.05). The number is the effect",
  "ratio relative to WT: odds ratio (larva-to-adult survival, hatching, Smurf),",
  "rate ratio (eggs per female) or time ratio (survival after a blood meal).",
  "* p < 0.05, ** p < 0.01, *** p < 0.001 (Dunnett-adjusted within line for",
  "eclosion, fecundity and hatching; Holm-adjusted within line for survival).",
  "† Ratio not estimable because one group had no events (no deaths or no",
  "Smurf+ females); significance from a log-rank or Fisher's exact test.",
  "‡ Not estimable: *cd*^g225^ heterozygous females died after blood",
  "feeding or did not feed. § Only 2 females. — Not assessed.",
  "WT groups: *cd*^g384^ eclosion and fecundity, non-transgenic siblings;",
  "*cd*^384R^ eclosion and fecundity, SDA-500; *cd*^g225^ fecundity, siblings",
  "pooled with SDA-500; survival after blood feeding, non-transgenic siblings",
  "of each line (*cd*^225R^ HET and trans-heterozygotes: *cd*^g225^ WT",
  "siblings). Transgenic mother = transgenic female × SDA-500 male;",
  "transgenic father = SDA-500 female × transgenic male.",
  "Observed counts, 95% CIs and exact p-values are given in Supplementary",
  "Table S_fitness."
)

# Pandoc markdown -> .docx (keeps italics and superscripts)
md_escape_cell <- \(x) str_replace_all(x, fixed("*"), "\\*")

table1_md <- c(
  "**Table 1. Fitness of *cardinal* knock-in and knock-out genotypes relative to wildtype.**",
  "",
  paste0(
    "| Group | Genotype | ",
    paste(paste(columns$trait, columns$sub, sep = ": "), collapse = " | "),
    " |"
  ),
  paste0("|", paste(rep("---", 2 + nrow(columns)), collapse = "|"), "|"),
  wide |>
    mutate(across(all_of(columns$col_key), md_escape_cell)) |>
    mutate(group = if_else(group == lag(group, default = ""), "", group)) |>
    glue::glue_data(
      "| {group} | {md} | ",
      "{eclosion_f} | {eclosion_m} | {eggs_mother} | {eggs_father} | ",
      "{hatch_mother} | {hatch_father} | {survival} | {smurf} |"
    ),
  "",
  table1_footnote
)
writeLines(table1_md, file.path(out_dir, "table1_fitness_summary.md"), useBytes = TRUE)

# HTML via gt
table1_gt <- wide |>
  select(group, html, all_of(columns$col_key)) |>
  gt(groupname_col = "group") |>
  fmt_markdown(columns = html) |>
  cols_label(
    html = "Genotype",
    eclosion_f = "Females",
    eclosion_m = "Males",
    eggs_mother = "Transgenic mother",
    eggs_father = "Transgenic father",
    hatch_mother = "Transgenic mother",
    hatch_father = "Transgenic father",
    survival = "Females",
    smurf = "Females"
  ) |>
  tab_spanner("Larva-to-adult survival", c(eclosion_f, eclosion_m)) |>
  tab_spanner("Fecundity (eggs per female)", c(eggs_mother, eggs_father)) |>
  tab_spanner("Fertility (hatching rate)", c(hatch_mother, hatch_father)) |>
  tab_spanner("Survival after blood meal", survival) |>
  tab_spanner("Midgut barrier loss (Smurf+)", smurf) |>
  cols_align("center", columns = all_of(columns$col_key)) |>
  tab_header(
    title = md("**Table 1. Fitness of *cardinal* knock-in and knock-out genotypes relative to wildtype**")
  ) |>
  tab_source_note(md(table1_footnote |> str_replace_all("\\^([^^]+)\\^", "<sup>\\1</sup>"))) |>
  opt_stylize(style = 1)

gtsave(table1_gt, file.path(out_dir, "table1_fitness_summary.html"))

# ── Supplementary table (detailed) ────────────────────────────────────────────

fmt_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    .default = formatC(p, format = "f", digits = 3)
  )
}
fmt_ratio_ci <- function(r, l, u) {
  f <- \(x) format(signif(x, 3), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  if_else(is.na(r), "not estimable", paste0(f(r), " [", f(l), "–", f(u), "]"))
}

table_s <- fitness_results |>
  mutate(
    row_key = factor(row_key, levels = rows$row_key),
    col_key = factor(col_key, levels = columns$col_key)
  ) |>
  arrange(row_key, col_key) |>
  transmute(
    group,
    genotype_md = md,
    genotype_html = html,
    trait,
    comparison = sub,
    observed = coalesce(observed, ""),
    effect = ratio_type,
    estimate_ci = if_else(
      is.na(p) & is.na(ratio),
      "not estimable",
      fmt_ratio_ci(ratio, lcl, ucl)
    ),
    p = fmt_p(p),
    test = coalesce(test, ""),
    note = coalesce(note, "")
  )

write_csv(
  table_s |> select(-genotype_html) |> rename(genotype = genotype_md),
  file.path(out_dir, "tableS_fitness_vs_wt.csv")
)

table_s_md <- c(
  "**Supplementary Table S_fitness. Fitness traits of each genotype compared with the wildtype group of the same line: observed values, effect ratio [95% CI] and p-value.**",
  "",
  "| Genotype | Trait | Comparison | Observed | Effect | Estimate [95% CI] | p | Test | Note |",
  "|---|---|---|---|---|---|---|---|---|",
  table_s |>
    mutate(note = str_replace_all(note, "cd\\[(\\w+)\\]", "*cd*^\\1^")) |>
    glue::glue_data(
      "| {genotype_md} | {trait} | {comparison} | {observed} | {effect} | ",
      "{estimate_ci} | {p} | {test} | {note} |"
    ),
  "",
  paste(
    "Observed values are raw counts or means from the data; estimates are from",
    "the fitted models (fecundity: zero-inflated negative binomial GLMM;",
    "hatching: binomial GLMM; larva-to-adult survival and Smurf: binomial GLM;",
    "survival: Weibull accelerated failure time model). Fecundity estimates are",
    "conditional means among non-structural-zero females, so they exceed the",
    "observed means where many females laid no eggs (proportion laying shown).",
    "Larva-to-adult counts are adults of each sex from 2 trays of 200",
    "mixed-sex larvae, so ~50% is the maximum under a 1:1 sex ratio.",
    "p-values: Dunnett-adjusted within line (eclosion, fecundity, hatching),",
    "Holm-adjusted within line (survival), or as stated; 95% CIs unadjusted."
  )
)
writeLines(table_s_md, file.path(out_dir, "tableS_fitness_vs_wt.md"), useBytes = TRUE)

table_s_gt <- table_s |>
  select(-genotype_md) |>
  gt(groupname_col = "group") |>
  fmt_markdown(columns = genotype_html) |>
  cols_label(
    genotype_html = "Genotype",
    trait = "Trait",
    comparison = "Comparison",
    observed = "Observed",
    effect = "Effect",
    estimate_ci = "Estimate [95% CI]",
    p = "p",
    test = "Test",
    note = "Note"
  ) |>
  tab_header(
    title = "Supplementary Table S_fitness. Fitness traits compared with wildtype of the same line"
  ) |>
  opt_stylize(style = 1)

gtsave(table_s_gt, file.path(out_dir, "tableS_fitness_vs_wt.html"))

# Convert the markdown tables to Word with Pandoc (bundled with Quarto)
for (f in c("table1_fitness_summary", "tableS_fitness_vs_wt")) {
  md_file <- file.path(out_dir, paste0(f, ".md"))
  status <- system2(
    "quarto",
    c("pandoc", shQuote(md_file), "-f", "markdown", "-o", shQuote(file.path(out_dir, paste0(f, ".docx"))))
  )
  if (status == 0) file.remove(md_file)
}

table1_gt
