source(here::here("scripts", "packages.R"))
source(here::here("scripts", "functions", "blood_feeding_functions.R"))


# ── Survival after blood feeding ───────────────────────────────────────────────

path <- here::here("data", "SurvivalAfterBFData.xlsx")

data_long <- path |>
  excel_sheets() |>
  set_names() |>
  map(\(s) read_excel(path = path, sheet = s)) |>
  list_rbind(names_to = "line") |>
  pivot_longer(cols = WT:Hom, names_to = "treatment", values_to = "event") |>
  drop_na(event) |>
  unite("line_treatment", c(line, treatment), remove = FALSE) |>
  mutate(line_treatment = fct_relevel(line_treatment, "D251_Hom"))

# Biological plausibility checks ----
# Hard checks: survival times must be positive and event indicator must be
# binary; verified before Weibull modelling (produces Figure 2)
data_long <- data_long |>
  verify(Hours > 0) |>             # survival time must be positive
  verify(event %in% c(0, 1))       # event must be binary censoring indicator

# Observed deaths per group (reported alongside model estimates)
bf_observed <- observed_deaths(data_long, line_treatment)

# Groups with no deaths (cd[g338-384] WT siblings, 0/47) cannot be modelled
# (complete separation) and, if left in, corrupt flexsurv's standard errors
# for every other group. Fit without them; report them as observed counts and
# show them in Figure 2 by their Kaplan-Meier curve only.
data_bf_model <- drop_zero_event_groups(data_long, line_treatment)
attr(data_bf_model, "zero_event_groups")

weibull_bf <- fit_weibull(data_bf_model, Surv(Hours, event) ~ line_treatment)
weibull_bf$summary

# Each Het/Hom group versus its own line's wildtype siblings (time ratios
# from the same Weibull model; Holm-adjusted within line)
bf_vs_wt <- tibble(level = levels(data_long$line_treatment)) |>
  filter(!str_detect(level, "_WT$")) |>
  mutate(
    ref = str_replace(level, "_(Het|Hom)$", "_WT"),
    family = str_remove(level, "_(Het|Hom)$")
  ) |>
  contrast_levels(model = weibull_bf$model, factor_name = "line_treatment") |>
  left_join(bf_observed, by = c("level" = "line_treatment")) |>
  left_join(
    bf_observed |> rename(n_ref = n, deaths_ref = deaths),
    by = c("ref" = "line_treatment")
  )
bf_vs_wt

surv_df <- pred_weibull(
  weibull_bf$model,
  data_bf_model,
  t_max = max(data_long$Hours)
) |>
  separate_wider_delim(
    line_treatment,
    delim = "_",
    names = c("line", "treatment"),
    too_many = "merge"
  ) |>
  mutate(lcl = if_else(est > 0.999, 1, lcl)) |> 
  filter(!(line =="1759" & treatment == "Hom")) # remove 1759 HOM due to low n


km_df <- tidy_km(data_long, Surv(Hours, event) ~ line_treatment) |>
 separate_wider_delim(strata, delim = "_", names = c("line", "treatment"), too_many = "merge") |> 
  filter(!(line =="1759" & treatment == "Hom")) # remove 1759 HOM due to low n





survival_after_blood <- ggplot(
  surv_df,
  aes(x = time, y = est, fill = factor(treatment))
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  geom_line(aes(colour = factor(treatment), linetype = lt_model)) +
  geom_step(
    data = km_df,
    aes(
      x = time,
      y = estimate,
      colour = factor(treatment),
      linetype = lt_km
    ),
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  scale_fill_manual(values = treatment_colours, name = "Treatment") +
  estimate_linetype_layers() +
  facet_wrap(~line, labeller = as_labeller(line_labels, label_parsed)) +
  labs(
    x = "Time (Hours)",
    y = "Survival Probability"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

survival_after_blood

# ── Trans-het survival ─────────────────────────────────────────────────────────

data_transhet <- read_excel(here::here("data", "TransHetSurvivalData.xlsx")) |>
  pivot_longer(
    cols = `2360Het`:`D251:2360WT`,
    names_to = "line_treatment",
    values_to = "event"
  ) |>
  drop_na(event) |>
  filter(Hours > 0)

# Column key: 2360Het = cd[g225] het; 2360WT = its wildtype siblings;
# D251:2360Het = cd[g225];cd[225R] trans-het; D251:2360WT = cd[225R] het
th_observed <- observed_deaths(data_transhet, line_treatment)

# cd[225R] hets had no deaths (0/46): fitted out for the same reason as in
# Figure 2 and reported as observed counts (Kaplan-Meier curve only in Fig 5)
data_th_model <- drop_zero_event_groups(data_transhet, line_treatment)
attr(data_th_model, "zero_event_groups")

weibull_th <- fit_weibull(data_th_model, Surv(Hours, event) ~ line_treatment)
weibull_th$summary

# All groups versus cd[g225] wildtype siblings, plus the key comparison of
# trans-hets versus cd[g225] hets (Holm-adjusted within each family)
th_contrasts <- tibble(
  level = c("2360Het", "D251:2360Het", "D251:2360WT", "D251:2360Het"),
  ref = c("2360WT", "2360WT", "2360WT", "2360Het"),
  family = c("vs_wt", "vs_wt", "vs_wt", "transhet_vs_ki")
) |>
  contrast_levels(model = weibull_th$model, factor_name = "line_treatment") |>
  left_join(th_observed, by = c("level" = "line_treatment")) |>
  left_join(
    th_observed |> rename(n_ref = n, deaths_ref = deaths),
    by = c("ref" = "line_treatment")
  )
th_contrasts

surv_df_trans <- pred_weibull(
  weibull_th$model,
  data_th_model,
  t_max = max(data_transhet$Hours) # trans-het follow-up ended at 108 h
) |>
  mutate(
    treatment = str_extract(line_treatment, "Het|WT"),
    line = str_remove(line_treatment, "Het|WT")
  ) |>
  mutate(lcl = if_else(est > 0.999, 1, lcl))

km_df_transhet <- tidy_km(data_transhet, Surv(Hours, event) ~ line_treatment) |>
  mutate(
    treatment = str_extract(strata, "Het|WT"),
    line = str_remove(strata, "Het|WT")
  )

# Panel specs: vary only legend title, labels, plot title, and y-axis
transhet_panels <- list(
  `2360` = list(
    legend_title = "<em>cd</em><sup>g225</sup>",
    legend_labels = c(
      "Het" = "<em>cd</em><sup>g225</sup> Het",
      "WT" = "<em>cd</em><sup>g225</sup> WT"
    ),
    plot_title = expression(italic(cd)^g225),
    colours = c("Het" = "#8BABD3", "WT" = "darkgrey"),
    show_y_label = TRUE,
    legend_order = 1
  ),
  `D251:2360` = list(
    legend_title = "<em>cd</em><sup>225R</sup>",
    legend_labels = c(
      "WT" = "<em>cd</em><sup>225R</sup> Het",
      "Het" = "<em>cd</em><sup>g225</sup>;<em>cd</em><sup>225R</sup> Het"
    ),
    plot_title = expression(
      italic(cd)^{
        225 * R
      }
    ),
    colours = c("Het" = "darkgreen", "WT" = "#FFA040"),
    show_y_label = FALSE,
    legend_order = 2
  )
)


transhet_plots <- imap(
  transhet_panels,
  ~ plot_transhet_panel(
    surv_df = surv_df_trans,
    km_df = km_df_transhet,
    line_val = .y,
    legend_title = NULL,
    legend_labels = .x$legend_labels,
    colours = .x$colours,
    plot_title = NULL,
    show_y_label = .x$show_y_label,
    legend_order = .x$legend_order
  )
)


# Collect legends so the shared Estimate key appears once, stacked to the
# right beneath the two genotype keys
transhet_survival <- transhet_plots$`2360` +
  transhet_plots$`D251:2360` +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.box = "vertical")
transhet_survival
