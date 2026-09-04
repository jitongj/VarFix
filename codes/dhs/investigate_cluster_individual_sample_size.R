# *******************
# date: 2026.07.28
# task: Investigate real DHS cluster-level individual sample sizes for Zambia
# author: Jitong Jiang
# ********************

library(dplyr)
library(ggplot2)
library(surveyPrev)
library(SUMMER)
library(rdhs)
library(here)


# **********************
# 1. Basic setting   -------
# **********************

indicator <- "CN_NUTS_C_WH2"
year <- 2018
country <- "Zambia"

dhsData <- getDHSdata(
  country = country,
  indicator = indicator,
  year = year
)

data <- getDHSindicator(
  dhsData,
  indicator = indicator
)

geo <- getDHSgeo(
  country = country,
  year = year
)

poly.adm1 <- readRDS(
  here::here("data", "poly.adm1.rds")
)

poly.adm2 <- readRDS(
  here::here("data", "poly.adm2.rds")
)

cluster.info <- clusterInfo(
  geo = geo,
  poly.adm1 = poly.adm1,
  poly.adm2 = poly.adm2,
  by.adm1 = "NAME_1",
  by.adm2 = "NAME_2"
)

data0 <- data %>%
  dplyr::filter(!is.na(value))
# **********************
# 2. Cluster-level individual sample size   -------
# **********************

cluster_sample_size <- data0 %>%
  dplyr::left_join(
    cluster.info$data,
    by = "cluster"
  ) %>%
  dplyr::filter(
    !is.na(value),
    !is.na(admin2.name.full)
  ) %>%
  dplyr::group_by(
    admin1.name,
    admin2.name.full,
    cluster,
    strata
  ) %>%
  dplyr::summarise(
    n_individual = dplyr::n(),
    .groups = "drop"
  )

# **********************
# 3. Fit Normal distribution   -------
# **********************

cluster_sample_size_trunc <- cluster_sample_size %>%
  dplyr::filter(
    n_individual >= 5,
    n_individual <= 30
  )

mu_n <- round(mean(cluster_sample_size_trunc$n_individual))
sd_n <- round(sd(cluster_sample_size_trunc$n_individual))

mu_n
sd_n




ggplot(
  cluster_sample_size,
  aes(x = n_individual)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 1,
    boundary = 0,
    color = "white"
  ) +
  stat_function(
    fun = dnorm,
    args = list(
      mean = mu_n,
      sd = sd_n
    ),
    linewidth = 1
  ) +
  labs(
    title = "Individual Sample Size per Cluster",
    subtitle = paste0(
      "Fitted Normal: mean = ", round(mu_n, 2),
      ", SD = ", round(sd_n, 2)
    ),
    x = "Number of individuals per cluster",
    y = "Density"
  ) +
  theme_bw()


# ******************************
# Cluster-level summary:-------
# sample size, sum of weights, and average individual weight
# ******************************

cluster_weight_summary <- data0 %>%
  dplyr::left_join(
    cluster.info$data,
    by = "cluster"
  ) %>%
  dplyr::filter(
    !is.na(value),
    !is.na(admin2.name.full),
    !is.na(weight)
  ) %>%
  dplyr::group_by(
    admin1.name,
    admin2.name.full,
    cluster,
    strata
  ) %>%
  dplyr::summarise(
    n_individual = dplyr::n(),
    sum_weight   = sum(weight, na.rm = TRUE),
    mean_weight  = mean(weight, na.rm = TRUE),
    .groups = "drop"
  )



## 1. sum of weights per cluster ----

summary_sum_weight <- cluster_weight_summary %>%
  dplyr::summarise(
    n_cluster = dplyr::n(),
    min       = min(sum_weight, na.rm = TRUE),
    q05       = quantile(sum_weight, 0.05, na.rm = TRUE),
    q10       = quantile(sum_weight, 0.10, na.rm = TRUE),
    q25       = quantile(sum_weight, 0.25, na.rm = TRUE),
    median    = median(sum_weight, na.rm = TRUE),
    mean      = mean(sum_weight, na.rm = TRUE),
    q75       = quantile(sum_weight, 0.75, na.rm = TRUE),
    q90       = quantile(sum_weight, 0.90, na.rm = TRUE),
    q95       = quantile(sum_weight, 0.95, na.rm = TRUE),
    max       = max(sum_weight, na.rm = TRUE)
  )

print(summary_sum_weight)



plot_sum_weight <- ggplot(
  cluster_weight_summary,
  aes(x = sum_weight)
) +
  geom_histogram(
    bins = 30,
    color = "white"
  ) +
  labs(
    title = "Distribution of Sum of Weights per Cluster",
    x = "Sum of individual weights in cluster",
    y = "Number of clusters"
  ) +
  theme_bw()

print(plot_sum_weight)



# 2. average individual weight per cluster -----
summary_mean_weight <- cluster_weight_summary %>%
  dplyr::summarise(
    n_cluster = dplyr::n(),
    min       = min(mean_weight, na.rm = TRUE),
    q05       = quantile(mean_weight, 0.05, na.rm = TRUE),
    q10       = quantile(mean_weight, 0.10, na.rm = TRUE),
    q25       = quantile(mean_weight, 0.25, na.rm = TRUE),
    median    = median(mean_weight, na.rm = TRUE),
    mean      = mean(mean_weight, na.rm = TRUE),
    q75       = quantile(mean_weight, 0.75, na.rm = TRUE),
    q90       = quantile(mean_weight, 0.90, na.rm = TRUE),
    q95       = quantile(mean_weight, 0.95, na.rm = TRUE),
    max       = max(mean_weight, na.rm = TRUE)
  )

print(summary_mean_weight)


plot_mean_weight <- ggplot(
  cluster_weight_summary,
  aes(x = mean_weight)
) +
  geom_histogram(
    bins = 30,
    color = "white"
  ) +
  labs(
    title = "Distribution of Average Individual Weight per Cluster",
    x = "Average individual weight in cluster",
    y = "Number of clusters"
  ) +
  theme_bw()

print(plot_mean_weight)



# 3. number of individuals per cluster
summary_cluster_size <- cluster_weight_summary %>%
  dplyr::summarise(
    n_cluster = dplyr::n(),
    min       = min(n_individual, na.rm = TRUE),
    q05       = quantile(n_individual, 0.05, na.rm = TRUE),
    q10       = quantile(n_individual, 0.10, na.rm = TRUE),
    q25       = quantile(n_individual, 0.25, na.rm = TRUE),
    median    = median(n_individual, na.rm = TRUE),
    mean      = mean(n_individual, na.rm = TRUE),
    q75       = quantile(n_individual, 0.75, na.rm = TRUE),
    q90       = quantile(n_individual, 0.90, na.rm = TRUE),
    q95       = quantile(n_individual, 0.95, na.rm = TRUE),
    max       = max(n_individual, na.rm = TRUE)
  )

print(summary_cluster_size)



plot_cluster_size <- ggplot(
  cluster_weight_summary,
  aes(x = n_individual)
) +
  geom_histogram(
    bins = 30,
    color = "white"
  ) +
  labs(
    title = "Distribution of Number of Individuals per Cluster",
    x = "Number of individuals in cluster",
    y = "Number of clusters"
  ) +
  theme_bw()

print(plot_cluster_size)
