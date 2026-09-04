# *******************
# date: 2026.01.20
# task: This script is for large sample data of Zambia, 
# author: Jitong Jiang
# ********************


library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(stringr)
library(surveyPrev)

# Load "large_sample50_nationalWeight.RData" and you can plot the figures
load(here::here("Data/large_sample50_nationalWeight.RData"))

# **********************
# 1. Basic setting   -------
# **********************
source(here::here("codes", "directEST_1030_national.R"))

indicator <- "CN_NUTS_C_WH2" 

year <- 2018
frame_year <-2010
country <- "Zambia"
country.abbrev = "ZMB"


dhsData <- getDHSdata(country = country, indicator = indicator, year = year)
data <- getDHSindicator(dhsData, indicator = indicator)
data0 <- data %>%
  filter(!is.na(value))


geo <- getDHSgeo(country = country, year = year)
poly.adm1 <- geodata::gadm(country=country.abbrev, level=1, path=tempdir())
poly.adm1 <- sf::st_as_sf(poly.adm1)
poly.adm2 <- geodata::gadm(country=country.abbrev, level=2, path=tempdir())
poly.adm2 <- sf::st_as_sf(poly.adm2) %>%
  mutate(admin2.name.full = paste0(NAME_1, "_", NAME_2))
# poly.adm2=poly.adm2[poly.adm2$ENGTYPE_2=="Local Authority",]


cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, 
                            by.adm1 = "NAME_1",by.adm2 = "NAME_2")
admin.info1 <- adminInfo(poly.adm = poly.adm1, admin = 1, by.adm = "NAME_1")
admin.info2 <- adminInfo(poly.adm = poly.adm2, admin = 2,
                         by.adm = "NAME_2", by.adm.upper = "NAME_1")


set.seed(2024)



# setting
m_true   <- 0.5
sigma1   <- 0.2
sigma2   <- 0.05
n_indiv  <- 50
n_iter   <- 1000
legal_tol <- 1e-12
conf_lvls <- c(0.50, 0.60, 0.70, 0.8, 0.9)

expit <- function(x) 1/(1+exp(-x))
logit <- function(p) log(p/(1-p))

# choose a Admin1
adm1 <- "Central"
adm2_tbl <- cluster.info[["data"]] %>%
  filter(admin1.name == adm1) %>%
  distinct(admin2.name, admin2.name.full)



# **********************
# 2. sampling function  -------
# **********************

gen_stratum_sample <- function(a1, a2, a2_full, st,
                               n_clusters = 50L,
                               n_indiv = n_indiv,
                               sigma1 = sigma1,
                               sigma2 = sigma2,
                               m_true = m_true) {
  # Generate n_clusters "selected clusters" for this strata, and assign a cluster random effect to each cluster.
  cluster_ids <- paste0(a2_full, "_", st, "_", seq_len(n_clusters))
  e_cluster   <- rnorm(n_clusters, mean = 0, sd = sigma1)
  
  # For each cluster, generate n_indiv individuals (Bernoulli) 
  purrr::pmap_dfr(
    list(cluster_ids, e_cluster),
    function(cl, e_c) {
      # Linear prediction: logit(m_true) + e_cluster + individual noise
      eta <- qlogis(m_true) + e_c + rnorm(n_indiv, 0, sigma2)
      tibble::tibble(
        admin1.name      = a1,
        admin2.name      = a2,
        admin2.name.full = a2_full,
        strata           = st,
        cluster          = cl,
        value            = rbinom(n_indiv, 1, plogis(eta)),
        weight           = 1L
      )
    }
  )
}

sample_once <- function() {
  purrr::map_dfr(seq_len(nrow(adm2_tbl)), function(i) {
    a2      <- adm2_tbl$admin2.name[i]
    a2_full <- adm2_tbl$admin2.name.full[i]
    
    # 50 clusters in each of the two layers
    urban_df <- gen_stratum_sample(
      a1 = adm1, a2 = a2, a2_full = a2_full, st = "urban",
      n_clusters = 50L, n_indiv = n_indiv,
      sigma1 = sigma1, sigma2 = sigma2, m_true = m_true
    )
    rural_df <- gen_stratum_sample(
      a1 = adm1, a2 = a2, a2_full = a2_full, st = "rural",
      n_clusters = 50L, n_indiv = n_indiv,
      sigma1 = sigma1, sigma2 = sigma2, m_true = m_true
    )
    
    dplyr::bind_rows(urban_df, rural_df)
  })
}


## 2.1df -----

compute_df_taylor_by_admin2 <- function(sampled){
  L_tab <- sampled %>%
    group_by(admin2.name.full, strata) %>%
    summarise(L_h = n_distinct(cluster), .groups="drop")
  L_tab %>%
    group_by(admin2.name.full) %>%
    summarise(K=sum(L_h[L_h>0], na.rm=TRUE),
              H=sum(L_h>0, na.rm=TRUE),
              df_unfixed=as.integer(K-H),
              df_fixed  =as.integer(K), .groups="drop")
}

## 2.2 one sim ----

process_one_draw <- function(sampled_df){
  df_tab <- compute_df_taylor_by_admin2(sampled_df)
  unfixed <- directEST_1030(sampled_df, NULL, admin=2, aggregation=FALSE,
                           var.fix=FALSE, all.fix=FALSE)$res.admin2 %>%
    mutate(type="unfixed")
  fixed_triggered_obj <- directEST_1030(sampled_df, NULL, admin=2, aggregation=FALSE,
                       var.fix=TRUE, all.fix=FALSE)
  trg <- fixed_triggered_obj$res.admin2 %>%
    mutate(type="fixed_triggered")
  fixed_area_names <- fixed_triggered_obj$fixed_areas
  allfx <- directEST_1030(sampled_df, NULL, admin=2, aggregation=FALSE,
                         var.fix=TRUE, all.fix=TRUE)$res.admin2 %>%
    mutate(type="fixed_all")

  area_sample_info <- sampled_df |>
    dplyr::group_by(admin2.name.full) |>
    dplyr::summarise(
      n_i = dplyr::n(),
      m_i = dplyr::n_distinct(cluster),
      .groups = "drop"
    )

  unfixed_benchmark <- unfixed %>%
    dplyr::left_join(area_sample_info, by = "admin2.name.full") %>%
    dplyr::mutate(
      admin2.name.full = trimws(admin2.name.full),
      legal = !(admin2.name.full %in% fixed_area_names)
    )

  fixed_input <- trg %>%
    dplyr::left_join(area_sample_info, by = "admin2.name.full") %>%
    dplyr::mutate(admin2.name.full = trimws(admin2.name.full)) %>%
    dplyr::left_join(
      unfixed_benchmark %>%
        dplyr::select(admin2.name.full, legal),
      by = "admin2.name.full"
    )

  gvf_data <- fixed_input %>%
    dplyr::mutate(
      log_var = log(direct.var),
      log_binom = log(direct.est * (1 - direct.est)),
      log_n = log(n_i),
      direct.var.fixed = direct.var
    ) %>%
    dplyr::filter(
      is.finite(log_var)
    )

  gvf_fit <- stats::lm(log_var ~ log_binom + log_n, data = gvf_data)

  gvf_data$gvf_log_var <- stats::predict(
    gvf_fit,
    newdata = gvf_data
  )
  gvf_data$gvf_var <- exp(gvf_data$gvf_log_var)
  gvf_data$gvf_logit_var <- gvf_data$gvf_var / (
    gvf_data$direct.est^2 *
      (1 - gvf_data$direct.est)^2
  )

  gvf_all <- gvf_data %>%
    dplyr::mutate(
      direct.var.gvf = gvf_var,
      direct.var = gvf_var,
      direct.logit.var = gvf_logit_var,
      direct.logit.prec = dplyr::if_else(
        is.na(direct.logit.var), NA_real_, 1 / direct.logit.var
      ),
      type = "GVF"
    )

  res <- bind_rows(unfixed, trg, allfx) %>%
    left_join(df_tab, by="admin2.name.full") %>%
    mutate(
      true_p = m_true,
      legal = !is.na(direct.var) & direct.var>legal_tol,
      df_use = if_else(type=="unfixed", df_unfixed, df_fixed)
    )
  
  all_unfixed <- res %>% filter(type=="unfixed") %>% mutate(method="all_unfixed")
  all_fixed   <- res %>% filter(type=="fixed_all") %>% mutate(method="all_fixed")
  fixed_illegal_unfixed_legal <- res %>% filter(type=="fixed_triggered") %>%
    mutate(method="fixed_illegal_unfixed_legal")
  unfixed_legal <- res %>% filter(type=="unfixed" & legal) %>%
    mutate(method="unfixed_legal")
  gvf_benchmark <- gvf_all %>%
    left_join(df_tab, by="admin2.name.full") %>%
    mutate(
      true_p = m_true,
      df_use = df_unfixed,
      method = "GVF"
    )
  bind_rows(all_unfixed, all_fixed, fixed_illegal_unfixed_legal, unfixed_legal, gvf_benchmark)
}

## 2.3 all sim ----
all_runs <- map_dfr(seq_len(n_iter), function(i){
  samp <- sample_once()
  out <- process_one_draw(samp)
  out$sim_id <- i
  out
})

# **********************
# 3. metrics function   -------
# **********************

summarise_metrics <- function(df_all, conf_levels=conf_lvls){
  one_level <- function(df, cl){
    alpha <- 1 - cl
    zlo <- qnorm((1 - cl)/2)
    zhi <- qnorm(1 - (1 - cl)/2)
    
    df %>%
     
      mutate(
        df_eff = case_when(
          method %in% c("all_unfixed", "fixed_illegal_unfixed_legal", "unfixed_legal", "GVF") ~ df_unfixed, # use K-H
          method %in% c("all_fixed") ~ df_fixed,                                                     # use K
          TRUE ~ df_use 
        ),
        se  = sqrt(pmax(direct.logit.var, 0)),
        tlo = if_else(!is.na(df_eff) & df_eff > 0, qt((1 - cl)/2, df_eff), zlo),
        thi = if_else(!is.na(df_eff) & df_eff > 0, qt(1 - (1 - cl)/2, df_eff), zhi),
        
        ci_lb_t = expit(direct.logit.est + tlo * se),
        ci_ub_t = expit(direct.logit.est + thi * se),
        ci_lb_z = expit(direct.logit.est + zlo * se),
        ci_ub_z = expit(direct.logit.est + zhi * se),
        
        covered_t = (true_p >= ci_lb_t & true_p <= ci_ub_t),
        covered_z = (true_p >= ci_lb_z & true_p <= ci_ub_z),
        
        width_t = ci_ub_t - ci_lb_t,
        width_z = ci_ub_z - ci_lb_z,
        
        score_t = width_t + (2/alpha) * (pmax(0, ci_lb_t - true_p) + pmax(0, true_p - ci_ub_t)),
        score_z = width_z + (2/alpha) * (pmax(0, ci_lb_z - true_p) + pmax(0, true_p - ci_ub_z)),
        
        sq_err = (direct.est - true_p)^2
      ) %>%
      group_by(admin2.name.full, method) %>%
      summarise(
        mse     = mean(sq_err, na.rm = TRUE),
        cov_t   = mean(covered_t, na.rm = TRUE),
        cov_z   = mean(covered_z, na.rm = TRUE),
        width_t = mean(width_t, na.rm = TRUE),
        width_z = mean(width_z, na.rm = TRUE),
        score_t = mean(score_t, na.rm = TRUE),
        score_z = mean(score_z, na.rm = TRUE),
        avg_var = mean(direct.var, na.rm = TRUE),
        emp_var = var(direct.est, na.rm = TRUE),
        var_disc= emp_var - avg_var,
        df_mean = mean(df_eff, na.rm = TRUE),   
        .groups = "drop"
      ) %>% mutate(conf_level = cl)
  }
  purrr::map_dfr(conf_levels, ~one_level(df_all, .x))
}

metrics <- summarise_metrics(all_runs)
print(head(metrics))

# **********************
# 4. Plotting  -------
# **********************

library(dplyr)
library(ggplot2)
library(forcats)
library(rlang)

## 3.1 Method Tags & Colors ----
.method3_levels <- c("all_unfixed","all_fixed","fixed_illegal_unfixed_legal","GVF")
.method3_labels <- c(
  all_unfixed = "All—Unfixed",
  all_fixed   = "All—Fixed",
  fixed_illegal_unfixed_legal = "Fixed(Illegal)+Unfixed(Legal)",
  GVF = "GVF"
)
.method3_colors <- c(
  "all_unfixed" = "#ff7f0e",
  "all_fixed"   = "#1f77b4",
  "fixed_illegal_unfixed_legal" = "#2ca02c",
  "GVF" = "#9467bd"
)
scale_method3_fill <- function() {
  cols <- setNames(.method3_colors, .method3_labels)
  scale_fill_manual(
    values = cols, 
    breaks = .method3_labels, 
    labels = .method3_labels,  
    name = "Method", 
    drop = FALSE
  )
}

## 3.2 Preprocessing: Retain comparison methods + sort by conf_lvls ----
prep_admin1_pool_faceted <- function(metrics, conf_lvls) {
  metrics %>%
    filter(method %in% .method3_levels) %>%
    mutate(
      method = factor(method, levels = .method3_levels, labels = .method3_labels),
      conf_level = factor(conf_level, levels = conf_lvls, labels = sprintf("%.2f", conf_lvls))
    )
}

## 3.3 General-purpose plotter: Admin1 summary + faceting (one cell per conf_level) ----
plot_admin1_box_faceted <- function(metrics, conf_lvls, y_col, title, ylab,
                                    clip01 = FALSE, add_target_line = FALSE) {
  pd <- prep_admin1_pool_faceted(metrics, conf_lvls)
  
  p <- ggplot(pd, aes(x = method, y = .data[[y_col]], fill = method)) +
    geom_boxplot(width = 0.6, alpha = 0.95, outlier.shape = NA) +
    facet_wrap(~conf_level, nrow = 1, scales = "fixed") +
    scale_method3_fill() +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 15, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
  
  if (clip01) {
    p <- p + coord_cartesian(ylim = c(0, 1))
  }
  
  if (add_target_line) {
    target_df <- tibble(
      conf_level = factor(sprintf("%.2f", conf_lvls), levels = levels(pd$conf_level)),
      yint = conf_lvls
    )
    p <- p + geom_hline(data = target_df, aes(yintercept = yint),
                        color = "red", linetype = "dashed", inherit.aes = FALSE)
  }
  
  p
}

## 3.4 8 images (Admin1 summary, three methods, faceted by conf_lvls) ----
# 1) Variance discrepancy
plot_admin1_var_disc_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "var_disc",
    title = "Variance Discrepancy",
    ylab  = "Empirical − Reported Variance",
    clip01 = FALSE, add_target_line = FALSE
  )
}

# 2) Coverage (Normal)
plot_admin1_cov_normal_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "cov_z",
    title = "Coverage",
    ylab  = "Empirical Coverage",
    clip01 = TRUE, add_target_line = TRUE
  )
}

# 3) Coverage (t)
plot_admin1_cov_t_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "cov_t",
    title = "Coverage (t) — Admin1",
    ylab  = "Empirical Coverage (t)",
    clip01 = TRUE, add_target_line = TRUE
  )
}

# 4) Degrees of Freedom (t)
plot_admin1_df_t_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "df_mean",
    title = "Degrees of Freedom (t) — Admin1",
    ylab  = "Mean df (per domain)",
    clip01 = FALSE, add_target_line = FALSE
  )
}

# 5) CI Width (Normal)
plot_admin1_width_normal_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "width_z",
    title = "CI Width",
    ylab  = "Average CI Width",
    clip01 = FALSE, add_target_line = FALSE
  )
}

# 6) CI Width (t)
plot_admin1_width_t_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "width_t",
    title = "CI Width (t)",
    ylab  = "Average CI Width (t)",
    clip01 = FALSE, add_target_line = FALSE
  )
}

# 7) Interval Score (Normal)
plot_admin1_is_normal_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "score_z",
    title = "Interval Score",
    ylab  = "Interval Score",
    clip01 = FALSE, add_target_line = FALSE
  )
}

# 8) Interval Score (t)
plot_admin1_is_t_faceted <- function(metrics, conf_lvls) {
  plot_admin1_box_faceted(
    metrics, conf_lvls,
    y_col = "score_t",
    title = "Interval Score (t)",
    ylab  = "Interval Score (t)",
    clip01 = FALSE, add_target_line = FALSE
  )
}


conf_lvls <- c(0.50, 0.60, 0.70, 0.8, 0.9) 
p_var   <- plot_admin1_var_disc_faceted(metrics, conf_lvls)
p_cov_z <- plot_admin1_cov_normal_faceted(metrics, conf_lvls)
p_cov_t <- plot_admin1_cov_t_faceted(metrics, conf_lvls)
p_df    <- plot_admin1_df_t_faceted(metrics, conf_lvls)
p_wz    <- plot_admin1_width_normal_faceted(metrics, conf_lvls)
p_wt    <- plot_admin1_width_t_faceted(metrics, conf_lvls)
p_sz    <- plot_admin1_is_normal_faceted(metrics, conf_lvls)
p_st    <- plot_admin1_is_t_faceted(metrics, conf_lvls)

print(p_cov_z);
print(p_wz);  print(p_sz)

# Save p_cov_z (Coverage - Normal)
ggsave(
  filename = here::here("figures", "sanity check",
                        "large_sample50_cover_nationalWeight.png"),
  plot = p_cov_z,
  width = 14, height = 6, dpi = 300
)

# Save p_wz (CI Width - Normal)
ggsave(
  filename = here::here("figures", "sanity check",
                        "large_sample50_ci_nationalWeight.png"),
  plot = p_wz,
  width = 14, height = 6, dpi = 300
)

# Save p_sz (Interval Score - Normal)
ggsave(
  filename = here::here("figures", "sanity check",
                        "large_sample50_score_nationalWeight.png"),
  plot = p_sz,
  width = 14, height = 6, dpi = 300
)

