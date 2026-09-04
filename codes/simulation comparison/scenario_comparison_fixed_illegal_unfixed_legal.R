# ********************
# Scenario comparison for Fixed(Illegal)+Unfixed(Legal)
# ********************

library(dplyr)
library(ggplot2)
library(here)

target_method <- "fixed_illegal_unfixed_legal"

scenario_levels <- c(
  "Default",
  "National Mean + 10% Weight",
  "Admin1 Mean + Admin1 Weight"
)

scenario_colors <- c(
  "Default" = "#2ca02c",
  "National Mean + 10% Weight" = "#1f77b4",
  "Admin1 Mean + Admin1 Weight" = "#9467bd"
)

scale_scenario_fill <- function() {
  scale_fill_manual(
    values = scenario_colors,
    breaks = scenario_levels,
    labels = scenario_levels,
    name = "Scenario",
    drop = FALSE
  )
}

default_metrics <- readRDS(here::here(
  "figures", "simulation", "direct_GVF_benchmarks",
  "simulation_direct_GVF_performance_summary.rds"
))

national_10w_metrics <- readRDS(here::here(
  "figures", "simulation_sensitivity", "nationalMean_10QWeight",
  "sensitivity10weight_simulation_direct_GVF_performance_summary.rds"
))

admin1_metrics <- readRDS(here::here(
  "figures", "simulation_sensitivity", "admin1Mean_admin1Weight",
  "sensitivityAdmin1MeanWeight_simulation_direct_GVF_performance_summary.rds"
))

all_summary_metrics <- bind_rows(
  default_metrics %>% mutate(scenario = "Default"),
  national_10w_metrics %>% mutate(scenario = "National Mean + 10% Weight"),
  admin1_metrics %>% mutate(scenario = "Admin1 Mean + Admin1 Weight")
) %>%
  filter(method == target_method) %>%
  mutate(scenario = factor(scenario, levels = scenario_levels))

plot_admin1_conf_08 <- function(metrics,
                                y_col,
                                title,
                                ylab,
                                add_target_line = FALSE) {
  conf_target <- 0.80

  pd <- metrics %>%
    filter(conf_level == conf_target)

  admin1_levels <- sort(unique(pd$admin1.name))

  pd <- pd %>%
    mutate(admin1.name = factor(admin1.name, levels = admin1_levels))

  p <- ggplot(
    pd,
    aes(x = admin1.name, y = .data[[y_col]], fill = scenario)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.7),
      width = 0.55,
      alpha = 0.95,
      outlier.shape = NA
    ) +
    scale_scenario_fill() +
    labs(title = title, x = "Admin1 Region", y = ylab) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  if (add_target_line) {
    p <- p + geom_hline(
      yintercept = conf_target,
      linetype = "dashed",
      color = "red"
    )
  }

  p
}

plot_admin1_bias_var_est <- function(metrics) {
  p <- plot_admin1_conf_08(
    metrics,
    y_col = "variance_discrepancy",
    title = "Bias of Variance Estimators",
    ylab = "Bias of Variance Estimators",
    add_target_line = FALSE
  )

  p + geom_hline(yintercept = 0, linetype = "dashed")
}

plot_admin1_cov_normal_08 <- function(metrics) {
  p <- plot_admin1_conf_08(
    metrics,
    y_col = "coverage_rate_z",
    title = "Coverage",
    ylab = "Empirical",
    add_target_line = TRUE
  )

  yvals <- metrics$coverage_rate_z
  y_min <- min(yvals, na.rm = TRUE)
  y_max <- max(yvals, na.rm = TRUE)

  p + coord_cartesian(ylim = c(y_min, y_max))
}

plot_admin1_width_normal_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "avg_width_z",
    title = "CI Width",
    ylab = "Average CI Width"
  )
}

plot_admin1_is_normal_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "avg_interval_score_z",
    title = "Interval Score",
    ylab = "Interval Score"
  )
}

p_bias_var_est <- plot_admin1_bias_var_est(all_summary_metrics)
p_cov_z <- plot_admin1_cov_normal_08(all_summary_metrics)
p_wz <- plot_admin1_width_normal_08(all_summary_metrics)
p_sz <- plot_admin1_is_normal_08(all_summary_metrics)

ggsave(
  here::here(
    "figures",
    "scenario_comparison_fixed_illegal_unfixed_legal_bias_variance.png"
  ),
  plot = p_bias_var_est,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  here::here(
    "figures",
    "scenario_comparison_fixed_illegal_unfixed_legal_coverage.png"
  ),
  plot = p_cov_z,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  here::here(
    "figures",
    "scenario_comparison_fixed_illegal_unfixed_legal_ci_width.png"
  ),
  plot = p_wz,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  here::here(
    "figures",
    "scenario_comparison_fixed_illegal_unfixed_legal_interval_score.png"
  ),
  plot = p_sz,
  width = 14,
  height = 6,
  dpi = 300
)
