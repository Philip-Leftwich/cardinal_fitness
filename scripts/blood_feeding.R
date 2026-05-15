source(here::here("scripts", "packages.R"))
source(here::here("scripts", "blood_feeding_functions.R"))


# ── Survival after blood feeding ───────────────────────────────────────────────

path <- "data/SurvivalAfterBFData.xlsx"

data_long <- path |>
  excel_sheets() |>
  set_names() |>
  map(\(s) read_excel(path = path, sheet = s)) |>
  list_rbind(names_to = "line") |>
  pivot_longer(cols = WT:Hom, names_to = "treatment", values_to = "event") |>
  drop_na(event) |>
  unite("line_treatment", c(line, treatment), remove = FALSE) |>
  mutate(line_treatment = fct_relevel(line_treatment, "D251_Hom"))

weibull_bf <- fit_weibull(data_long, Surv(Hours, event) ~ line_treatment)
weibull_bf$summary

surv_df <- pred_weibull(
  weibull_bf$model,
  data_long,
  t_max = max(data_long$Hours)
) |>
  separate(
    line_treatment,
    into = c("line", "treatment"),
    sep = "_",
    extra = "merge"
  ) |>
  mutate(lcl = if_else(est > 0.999, 1, lcl))


km_df <- tidy_km(data_long, Surv(Hours, event) ~ line_treatment) |>
  separate(strata, into = c("line", "treatment"), sep = "_", extra = "merge")

survival_after_blood <- ggplot(
  surv_df,
  aes(x = time, y = est, fill = factor(treatment))
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  geom_line(aes(colour = factor(treatment))) +
  geom_step(
    data = km_df,
    aes(x = time, y = estimate, colour = factor(treatment)),
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = treatment_colours, name = "Treatment") +
  scale_fill_manual(values = treatment_colours, name = "Treatment") +
  facet_wrap(~line, labeller = as_labeller(line_labels, label_parsed)) +
  labs(
    x = "Time (Hours)",
    y = "Survival Probability",
    caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

# ── Trans-het survival ─────────────────────────────────────────────────────────

data_transhet <- read_excel("data/TransHetSurvivalData.xlsx") |>
  pivot_longer(
    cols = `2360Het`:`D251:2360WT`,
    names_to = "line_treatment",
    values_to = "event"
  ) |>
  drop_na(event) |>
  filter(Hours > 0)

weibull_th <- fit_weibull(data_transhet, Surv(Hours, event) ~ line_treatment)
weibull_th$summary

surv_df_trans <- pred_weibull(
  weibull_th$model,
  data_transhet,
  t_max = max(data_long$Hours)
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
    show_y_label = TRUE
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
    show_y_label = FALSE
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
    show_y_label = .x$show_y_label
  )
)


transhet_plots$`D251:2360` <- transhet_plots$`D251:2360` +
  labs(caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier")

transhet_survival <- transhet_plots$`2360` + transhet_plots$`D251:2360`
transhet_survival
