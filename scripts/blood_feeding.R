source(here::here("scripts", "packages.R"))

# ── Shared lookups ─────────────────────────────────────────────────────────────

line_labels <- c(
  "1759" = "italic(cd)^g384",
  "2360B5" = "italic(cd)^g225",
  "2072" = "italic(cd)^g384_del",
  "2301" = "italic(cd)^{g338-384}",
  "D251" = "italic(cd)^{225*R}",
  "QA383P" = "italic(cd)^{384*R}"
)

treatment_colours <- c(
  "WT" = "#4D4D4D",
  "Het" = "#E69F00",
  "Hom" = "#56B4E9"
)

# ── Helper functions ───────────────────────────────────────────────────────────

# Fit Weibull model and return tidy HR summary
fit_weibull <- function(data, formula) {
  model <- flexsurvreg(formula, data = data, dist = "weibull")
  model$call$data <- substitute(data)
  model$call$formula <- substitute(formula)
  summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  list(model = model, summary = summary)
}

# Predict survival from a fitted flexsurvreg model
pred_weibull <- function(model, data, t_max, n_t = 108) {
  pred_data <- expand_grid(line_treatment = unique(data$line_treatment))
  summary(
    model,
    newdata = pred_data,
    type = "survival",
    t = seq(0, t_max, length.out = n_t),
    tidy = TRUE
  ) |>
    as_tibble()
}

# Tidy Kaplan-Meier fit, stripping strata prefix
tidy_km <- function(data, formula) {
  km_fit <- survfit(formula, data = data)
  strata_prefix <- paste0(all.vars(formula[[3]]), "=")
  broom::tidy(km_fit) |>
    mutate(
      strata = str_remove_all(strata, paste(strata_prefix, collapse = "|"))
    )
}

# Build one transhet panel: filter surv_df/km_df by line_val, apply custom scales
plot_transhet_panel <- function(
  surv_df,
  km_df,
  line_val,
  legend_title,
  legend_labels,
  plot_title,
  show_y_label = TRUE
) {
  p <- ggplot(
    filter(surv_df, line == line_val),
    aes(
      x = time,
      y = est,
      colour = treatment,
      linetype = treatment,
      fill = treatment
    )
  ) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_step(
      data = filter(km_df, line == line_val),
      aes(x = time, y = estimate, colour = treatment, linetype = treatment),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = c("Het" = "#E69F00", "WT" = "#4D4D4D"),
      name = legend_title,
      labels = legend_labels
    ) +
    scale_fill_manual(
      values = c("Het" = "#E69F00", "WT" = "#4D4D4D"),
      name = legend_title,
      labels = legend_labels
    ) +
    scale_linetype_manual(
      values = c("Het" = "solid", "WT" = "dashed"),
      name = legend_title,
      labels = legend_labels
    ) +
    ggtitle(plot_title) +
    transhet_theme

  if (!show_y_label) {
    p <- p + labs(y = NULL)
  }
  p
}

# Shared theme for transhet panels
transhet_theme <- list(
  scale_y_continuous(limits = c(0, 1)),
  labs(
    x = "Time (Hours)",
    y = "Survival Probability",
    caption = "Solid: Weibull model  |  Dashed: Kaplan-Meier"
  ),
  theme_bw(base_size = 12),
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    panel.grid.minor = element_blank()
  )
)

# ── Survival after blood feeding ───────────────────────────────────────────────

path <- "clean_data/SurvivalAfterBFData.xlsx"

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
  )

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

data_transhet <- read_excel("clean_data/TransHetSurvivalData.xlsx") |>
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
  )

km_df_transhet <- tidy_km(data_transhet, Surv(Hours, event) ~ line_treatment) |>
  mutate(
    treatment = str_extract(strata, "Het|WT"),
    line = str_remove(strata, "Het|WT")
  )

# Panel specs: vary only legend title, labels, plot title, and y-axis
transhet_panels <- list(
  `2360` = list(
    legend_title = "<em>cd</em><sup>g225</sup>",
    legend_labels = c("Het" = "Het", "WT" = "WT"),
    plot_title = expression(italic(cd)^g225),
    show_y_label = TRUE
  ),
  `D251:2360` = list(
    legend_title = "<em>cd</em><sup>225R</sup>",
    legend_labels = c("Het" = "Het<sup>KO</sup>", "WT" = "Trans-het"),
    plot_title = expression(
      italic(cd)^{
        225 * R
      }
    ),
    show_y_label = FALSE
  )
)

transhet_plots <- imap(
  transhet_panels,
  ~ plot_transhet_panel(
    surv_df = surv_df_trans,
    km_df = km_df_transhet,
    line_val = .y,
    legend_title = .x$legend_title,
    legend_labels = .x$legend_labels,
    plot_title = .x$plot_title,
    show_y_label = .x$show_y_label
  )
)

transhet_survival <- transhet_plots$`2360` + transhet_plots$`D251:2360`
