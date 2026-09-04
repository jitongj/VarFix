# *******************
# date: 2026.08.09
# task: Zambia admin2 real-data direct, FH, GVF, GVF-FH, and VSALM analysis.
# note: This is real-data analysis, not simulation.
# ********************

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)
library(surveyPrev)
library(SUMMER)
library(here)

source(here::here("codes", "directEST_1030_national.R"))
source(here::here("codes", "fhModel_1030.R"))

vsalm_local_path <- here::here("VSALM-paper-main", "VSALM-main")
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(vsalm_local_path, quiet = TRUE)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(vsalm_local_path, quiet = TRUE)
} else {
  stop("Install devtools or pkgload to load the local modified VSALM package.")
}

## **************************
## Basic settings and Zambia real DHS data ----------
## **************************

country <- "Zambia"
country.abbrev <- "ZMB"
indicator <- "CN_NUTS_C_WH2"
year <- 2018
ci_level <- 0.95

dhsData <- getDHSdata(
  country = country,
  indicator = indicator,
  year = year
)

data <- getDHSindicator(
  dhsData,
  indicator = indicator
)

data0 <- data %>%
  dplyr::filter(!is.na(value))

geo <- getDHSgeo(
  country = country,
  year = year
)

poly.adm1 <- readRDS(
  here::here("data", "poly.adm1.rds")
)

poly.adm2 <- readRDS(
  here::here("data", "poly.adm2.rds")
) %>%
  dplyr::mutate(
    admin2.name.full = paste0(NAME_1, "_", NAME_2)
  )

cluster.info <- clusterInfo(
  geo = geo,
  poly.adm1 = poly.adm1,
  poly.adm2 = poly.adm2,
  by.adm1 = "NAME_1",
  by.adm2 = "NAME_2"
)

admin.info2 <- adminInfo(
  poly.adm = poly.adm2,
  admin = 2,
  by.adm = "NAME_2",
  by.adm.upper = "NAME_1"
)

dat <- data0 %>%
  dplyr::left_join(
    cluster.info$data,
    by = "cluster"
  ) %>%
  dplyr::filter(
    !is.na(admin2.name.full)
  ) %>%
  dplyr::select(
    value,
    admin1.name,
    admin2.name,
    admin2.name.full,
    strata,
    cluster,
    householdID,
    weight
  )

## **************************
## Admin2 direct estimates ----------
## **************************

unfixed_admin2 <- directEST_1030(
  data = dat,
  NULL,
  admin = 2,
  CI = ci_level,
  aggregation = FALSE,
  var.fix = FALSE,
  all.fix = FALSE
)

unfixed <- unfixed_admin2$res.admin2 %>%
  dplyr::mutate(
    method = "direct_unfixed"
  )

unfixed_admin1 <- directEST_1030(
  data = dat,
  NULL,
  admin = 1,
  CI = ci_level
)

fixed_admin2 <- directEST_1030(
  data = dat,
  NULL,
  admin = 2,
  CI = ci_level,
  aggregation = FALSE,
  var.fix = TRUE,
  all.fix = FALSE
)

fixed <- fixed_admin2$res.admin2 %>%
  dplyr::mutate(
    method = "direct_fixed"
  )

fixed_area_names <- fixed_admin2$fixed_areas

fixed <- fixed %>%
  dplyr::mutate(
    original_status = dplyr::if_else(
      admin2.name.full %in% fixed_area_names,
      "Originally illegal",
      "Originally legal"
    )
  )

## **************************
## Admin2 nested FH estimates from interval plot framework ----------
## **************************

bad_clusters <- cluster.info$data %>%
  dplyr::filter(admin2.name.full %in% fixed_area_names) %>%
  dplyr::pull(cluster)

fh_unfixed_fit <- fhModel_1030(
  data = data0 %>%
    dplyr::filter(!cluster %in% bad_clusters),
  cluster.info = cluster.info,
  admin.info = admin.info2,
  admin = 2,
  CI = ci_level,
  model = "bym2",
  aggregation = FALSE,
  var.fix = FALSE,
  nested = TRUE
)

fh_fixed_fit <- fhModel_1030(
  data = data0,
  cluster.info = cluster.info,
  admin.info = admin.info2,
  admin = 2,
  CI = ci_level,
  model = "bym2",
  aggregation = FALSE,
  var.fix = TRUE,
  nested = TRUE
)

fh_unfixed <- fh_unfixed_fit$res.admin2 %>%
  dplyr::mutate(method = "FH")

fh_fixed <- fh_fixed_fit$res.admin2 %>%
  dplyr::mutate(method = "Fixed_FH")

## **************************
## Validate fixed direct estimates and variances ----------
## **************************

point_check <- all(
  is.finite(fixed$direct.est) &
    fixed$direct.est > 0 &
    fixed$direct.est < 1
)

variance_check <- all(
  is.finite(fixed$direct.var) &
    fixed$direct.var > 0 &
    is.finite(fixed$direct.logit.var) &
    fixed$direct.logit.var > 0
)

if (!point_check) {
  warning(
    paste0(
      "Some triggered-fixed direct estimates are non-finite ",
      "or outside (0, 1). GVF and logit VSALM results may be invalid."
    )
  )
}

if (!variance_check) {
  warning(
    paste0(
      "Some triggered-fixed direct variances are non-finite ",
      "or non-positive. GVF and VSALM results may be invalid."
    )
  )
}

## **************************
## Common admin2-level input table ----------
## **************************

area_sample_info <- dat %>%
  dplyr::group_by(
    admin2.name.full,
    admin1.name
  ) %>%
  dplyr::summarise(
    n_i = dplyr::n(),
    m_i = dplyr::n_distinct(cluster),
    H_i = dplyr::n_distinct(strata),
    .groups = "drop"
  )

fixed_input <- fixed %>%
  dplyr::select(
    admin2.name.full,
    admin1.name,
    admin2.name,
    direct.est,
    direct.var,
    direct.logit.est,
    direct.logit.var,
    direct.lower,
    direct.upper,
    original_status
  ) %>%
  dplyr::left_join(
    area_sample_info,
    by = c(
      "admin2.name.full",
      "admin1.name"
    )
  ) %>%
  dplyr::mutate(
    df = if_else(
      original_status == "Originally illegal",
      m_i,
      pmax(m_i - H_i, 1L)
    ),
    fixed_lower = direct.lower,
    fixed_upper = direct.upper,
    fixed_ci_width = direct.upper -
      direct.lower
  )

df_check <- all(
  is.finite(fixed_input$df) &
    fixed_input$df > 0
)

if (!df_check) {
  warning(
    "Some admin2 areas have non-positive VSALM degrees of freedom."
  )
}

df_problem_areas <- fixed_input %>%
  dplyr::filter(
    !is.finite(df) |
      df <= 0
  ) %>%
  dplyr::select(
    admin2.name.full,
    admin1.name,
    m_i,
    H_i,
    df
  )

## **************************
## GVF  ----------
## **************************

gvf_data <- fixed_input %>%
  dplyr::mutate(
    log_var = log(direct.var),
    log_binom = log(
      direct.est * (1 - direct.est)
    ),
    log_n = log(n_i)
  )

gvf_fit <- stats::lm(
  log_var ~ log_binom + log_n,
  data = gvf_data
)

gvf_data$gvf_log_var <- stats::predict(
  gvf_fit,
  newdata = gvf_data
)

gvf_data$gvf_var <- exp(
  gvf_data$gvf_log_var
)

gvf_data$gvf_logit_var <- gvf_data$gvf_var / (
  gvf_data$direct.est^2 *
    (1 - gvf_data$direct.est)^2
)

z_lower <- stats::qnorm((1 - ci_level) / 2)
z_upper <- stats::qnorm(1 - (1 - ci_level) / 2)

gvf_results <- gvf_data %>%
  dplyr::transmute(
    admin2.name.full,
    admin1.name,
    original_status,
    n_i,
    m_i,
    H_i,
    df,
    
    # GVF does not change the point estimate
    direct.est,
    
    fixed_var = direct.var,
    gvf_var,
    gvf_log_var,
    
    fixed_logit_var = direct.logit.var,
    gvf_logit_var,
    lower = SUMMER::expit(
      direct.logit.est + z_lower * sqrt(gvf_logit_var)
    ),
    upper = SUMMER::expit(
      direct.logit.est + z_upper * sqrt(gvf_logit_var)
    ),
    interval_width = upper - lower,
    
    method = "GVF_fixed"
  )

## **************************
## Admin2 adjacency matrix for VSALM BYM2 ----------
## **************************

vsalm_input <- admin.info2$data %>%
  dplyr::select(
    admin1.name,
    admin2.name,
    admin2.name.full
  ) %>%
  dplyr::left_join(
    fixed_input %>%
      dplyr::select(
        admin2.name.full,
        direct.est,
        direct.var,
        direct.logit.est,
        direct.logit.var,
        direct.lower,
        direct.upper,
        original_status,
        n_i,
        m_i,
        H_i,
        df,
        fixed_lower,
        fixed_upper,
        fixed_ci_width
      ),
    by = "admin2.name.full"
  ) %>%
  dplyr::mutate(
    original_status = dplyr::if_else(
      is.na(original_status),
      "Missing",
      original_status
    )
  )

poly_match <- match(
  vsalm_input$admin2.name.full,
  poly.adm2$admin2.name.full
)

if (any(is.na(poly_match))) {
  stop(
    paste0(
      "Some VSALM admin2 areas are missing from poly.adm2; ",
      "cannot build aligned admin2 adjacency matrix."
    )
  )
}

poly.adm2_vsalm <- poly.adm2 %>%
  dplyr::slice(poly_match)

admin2_adj_nb <- spdep::poly2nb(
  poly.adm2_vsalm,
  row.names = vsalm_input$admin2.name.full
)

admin2_adj_mat <- spdep::nb2mat(
  admin2_adj_nb,
  style = "B",
  zero.policy = TRUE
)

rownames(admin2_adj_mat) <- vsalm_input$admin2.name.full
colnames(admin2_adj_mat) <- vsalm_input$admin2.name.full

## **************************
## VSALM BYM2 spatial joint logit model ----------
## **************************

vsalm_var_tol <- 0

fit_with_warnings <- function(expr) {
  
  warnings_seen <- character()
  
  value <- withCallingHandlers(
    tryCatch(
      expr,
      error = function(e) e
    ),
    warning = function(w) {
      warnings_seen <<- c(
        warnings_seen,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )
  
  list(
    value = value,
    warnings = unique(warnings_seen),
    error = if (inherits(value, "error")) {
      conditionMessage(value)
    } else {
      NA_character_
    }
  )
}

## GVF-FH model

gvf_fh_direct_logit <- gvf_data %>%
  dplyr::transmute(
    region = admin2.name.full,
    HT.logit.est = direct.logit.est,
    HT.logit.var = dplyr::if_else(
      is.finite(gvf_logit_var) & gvf_logit_var > 0,
      gvf_logit_var,
      NA_real_
    )
  )

gvf_fh_missing_regions <- setdiff(
  admin.info2$data$admin2.name.full,
  gvf_fh_direct_logit$region
)

if (length(gvf_fh_missing_regions) > 0) {
  gvf_fh_direct_logit <- dplyr::bind_rows(
    gvf_fh_direct_logit,
    tibble::tibble(
      region = gvf_fh_missing_regions,
      HT.logit.est = NA_real_,
      HT.logit.var = NA_real_
    )
  )
}

gvf_fh_X <- admin.info2$data %>%
  dplyr::select(
    region = admin2.name.full,
    admin1.name
  )

gvf_fh_fit_obj <- fit_with_warnings(
  SUMMER::smoothSurvey(
    data = NULL,
    direct.est = gvf_fh_direct_logit,
    Amat = admin.info2$mat,
    regionVar = "region",
    responseVar = "HT.logit.est",
    direct.est.var = "HT.logit.var",
    response.type = "gaussian",
    CI = ci_level,
    smooth = TRUE,
    save.draws = TRUE,
    X = gvf_fh_X
  )
)

gvf_fh_fit <- gvf_fh_fit_obj$value

gvf_fh_results <- if (inherits(gvf_fh_fit, "error")) {
  tibble::tibble()
} else {
  gvf_fh_fit$smooth$logit.mean <- gvf_fh_fit$smooth$mean
  gvf_fh_fit$smooth$logit.var <- gvf_fh_fit$smooth$var
  gvf_fh_fit$smooth$logit.median <- gvf_fh_fit$smooth$median
  gvf_fh_fit$smooth$logit.lower <- gvf_fh_fit$smooth$lower
  gvf_fh_fit$smooth$logit.upper <- gvf_fh_fit$smooth$upper
  
  gvf_fh_prob_draws <- apply(
    gvf_fh_fit$draws.est[, -c(1, 2)],
    2,
    SUMMER::expit
  )
  
  gvf_fh_fit$smooth$mean <- apply(gvf_fh_prob_draws, 1, mean)
  gvf_fh_fit$smooth$var <- apply(gvf_fh_prob_draws, 1, var)
  gvf_fh_fit$smooth$median <- apply(gvf_fh_prob_draws, 1, median)
  gvf_fh_fit$smooth$lower <- apply(
    gvf_fh_prob_draws,
    1,
    stats::quantile,
    probs = (1 - ci_level) / 2
  )
  gvf_fh_fit$smooth$upper <- apply(
    gvf_fh_prob_draws,
    1,
    stats::quantile,
    probs = 1 - (1 - ci_level) / 2
  )
  
  gvf_fh_fit$smooth %>%
    dplyr::rename(admin2.name.full = region) %>%
    dplyr::left_join(
      admin.info2$data %>%
        dplyr::select(
          admin1.name,
          admin2.name,
          admin2.name.full
        ),
      by = "admin2.name.full"
    ) %>%
    dplyr::transmute(
      admin2.name.full,
      admin1.name,
      admin2.name,
      estimate = mean,
      variance = var,
      lower,
      upper,
      interval_width = upper - lower,
      method = "GVF_FH"
    )
}

## BYM2 model

inla.qinv <- INLA::inla.qinv

vsalm_bym2_fit_obj <- fit_with_warnings(
  VSALM::spatialJointSmoothLogit(
    Yhat = vsalm_input$direct.est,
    Vhat = vsalm_input$direct.var,
    domain = vsalm_input$admin2.name.full,
    na = as.integer(vsalm_input$n_i),
    df = as.integer(vsalm_input$df),
    adj_mat = admin2_adj_mat,
    var_tol = vsalm_var_tol,
    detailed_output = TRUE
  )
)

vsalm_bym2_fit <- vsalm_bym2_fit_obj$value

## **************************
## Extract VSALM estimates and credible intervals ----------
## **************************

vsalm_bym2_results <- if (inherits(vsalm_bym2_fit, "error")) {
  tibble::tibble()
} else {
  vsalm_bym2_fit$est %>%
    dplyr::transmute(
      admin2.name.full = domain,
      posterior_mean = mean,
      posterior_variance = var,
      lower,
      upper,
      posterior_ci_width = upper - lower,
      method = "VSALM-BYM2"
    )
}

## **************************
## Combine VSALM results ----------
## **************************

vsalm_results <- dplyr::bind_rows(
  vsalm_bym2_results
) %>%
  dplyr::left_join(
    vsalm_input %>%
      dplyr::select(
        admin2.name.full,
        admin1.name,
        original_status,
        n_i,
        m_i,
        H_i,
        df,
        direct.est,
        direct.var,
        fixed_lower,
        fixed_upper,
        fixed_ci_width
      ),
    by = "admin2.name.full"
  )

## **************************
## Results and diagnostics ----------
## **************************

direct_results <- dplyr::bind_rows(
  unfixed,
  fixed
)

common_area_info <- admin.info2$data %>%
  dplyr::transmute(
    area = admin2.name.full,
    admin1.name,
    admin2.name
  ) %>%
  dplyr::left_join(
    fixed_input %>%
      dplyr::transmute(
        area = admin2.name.full,
        original_status,
        n_i,
        m_i,
        H_i,
        df
      ),
    by = "area"
  ) %>%
  dplyr::mutate(
    original_status = dplyr::if_else(
      is.na(original_status),
      "Missing",
      original_status
    )
  )

all_method_results <- dplyr::bind_rows(
  unfixed %>%
    dplyr::transmute(
      area = admin2.name.full,
      method = "Direct",
      estimate = direct.est,
      variance = direct.var,
      lower = direct.lower,
      upper = direct.upper
    ),
  fixed %>%
    dplyr::transmute(
      area = admin2.name.full,
      method = "Fix Direct",
      estimate = direct.est,
      variance = direct.var,
      lower = direct.lower,
      upper = direct.upper
    ),
  fh_unfixed %>%
    dplyr::transmute(
      area = admin2.name.full,
      method = "FH",
      estimate = mean,
      variance = var,
      lower,
      upper
    ),
  fh_fixed %>%
    dplyr::transmute(
      area = admin2.name.full,
      method = "Fixed FH",
      estimate = mean,
      variance = var,
      lower,
      upper
    ),
  gvf_results %>%
    dplyr::transmute(
      area = admin2.name.full,
      method = "GVF",
      estimate = direct.est,
      variance = gvf_var,
      lower,
      upper
    ),
  gvf_fh_results %>%
    dplyr::transmute(
      area = admin2.name.full,
      method = "GVF-FH",
      estimate,
      variance,
      lower,
      upper
    ),
  vsalm_results %>%
    dplyr::transmute(
      area = admin2.name.full,
      method,
      estimate = posterior_mean,
      variance = posterior_variance,
      lower,
      upper
    )
) %>%
  dplyr::left_join(
    common_area_info,
    by = "area"
  ) %>%
  dplyr::mutate(
    interval_width = upper - lower,
    coverage = NA_real_,
    method = factor(
      method,
      levels = c(
        "Direct",
        "Fix Direct",
        "FH",
        "Fixed FH",
        "GVF",
        "GVF-FH",
        "VSALM-BYM2"
      )
    )
  )

method_results <- list(
  direct_unfixed = unfixed,
  direct_fixed = fixed,
  fh_unfixed = fh_unfixed,
  fh_fixed = fh_fixed,
  gvf_fixed = gvf_results,
  gvf_fh = gvf_fh_results,
  vsalm_bym2 = vsalm_bym2_results,
  vsalm_combined = vsalm_results,
  all_methods = all_method_results
)

analysis_diagnostics <- list(
  point_check = point_check,
  variance_check = variance_check,
  df_check = df_check,
  df_problem_areas = df_problem_areas,
  fixed_area_names = fixed_area_names,
  vsalm_var_tol = vsalm_var_tol,
  gvf_training_area_n = nrow(gvf_data),
  vsalm_bym2_area_n = nrow(vsalm_bym2_results),
  admin2_adj_mat_dim = dim(admin2_adj_mat),
  gvf_fh_warnings = gvf_fh_fit_obj$warnings,
  gvf_fh_error = gvf_fh_fit_obj$error,
  vsalm_bym2_warnings = vsalm_bym2_fit_obj$warnings,
  vsalm_bym2_error = vsalm_bym2_fit_obj$error
)

## **************************
## GVF plot data ----------
## **************************

gvf_point_scatter_data <- gvf_results %>%
  dplyr::transmute(
    admin2.name.full,
    admin1.name,
    original_status,
    n_i,
    fixed_estimate = direct.est,
    gvf_estimate = direct.est
  )

gvf_variance_scatter_data <- gvf_results %>%
  dplyr::transmute(
    admin2.name.full,
    admin1.name,
    original_status,
    n_i,
    fixed_logit_variance = fixed_logit_var,
    gvf_logit_variance = gvf_logit_var
  )

## **************************
## VSALM plot data ----------
## **************************

vsalm_method_levels <- c(
  "VSALM-BYM2"
)

status_levels <- c(
  "Originally legal",
  "Originally illegal",
  "Missing"
)

vsalm_point_scatter_data <- vsalm_results %>%
  dplyr::transmute(
    admin2.name.full,
    admin1.name,
    original_status,
    n_i,
    method,
    fixed_estimate = direct.est,
    vsalm_estimate = posterior_mean
  ) %>%
  dplyr::mutate(
    method = factor(
      method,
      levels = vsalm_method_levels
    ),
    original_status = factor(
      original_status,
      levels = status_levels
    )
  )

vsalm_ci_width_scatter_data <- vsalm_results %>%
  dplyr::transmute(
    admin2.name.full,
    admin1.name,
    original_status,
    n_i,
    method,
    fixed_ci_width = fixed_ci_width,
    vsalm_ci_width = posterior_ci_width
  ) %>%
  dplyr::mutate(
    method = factor(
      method,
      levels = vsalm_method_levels
    ),
    original_status = factor(
      original_status,
      levels = status_levels
    )
  )

method_results$gvf_point_comparison <- gvf_point_scatter_data
method_results$gvf_logit_variance_comparison <- gvf_variance_scatter_data
method_results$vsalm_point_comparison <- vsalm_point_scatter_data
method_results$vsalm_ci_width_comparison <- vsalm_ci_width_scatter_data

## **************************
## Combined interval plot and comparison data ----------
## **************************

combined_data <- all_method_results %>%
  dplyr::transmute(
    admin1.name,
    admin2.name.full = area,
    estimate,
    lower,
    upper,
    model = dplyr::if_else(
      as.character(method) == "VSALM-BYM2",
      "Joint",
      as.character(method)
    ),
    area_type = dplyr::case_when(
      area %in% fixed_area_names ~ "illegal variance",
      original_status == "Missing" ~ "missing",
      TRUE ~ "legal variance"
    )
  ) %>%
  dplyr::mutate(
    area_type = factor(
      area_type,
      levels = c("legal variance", "illegal variance", "missing")
    ),
    model = factor(
      model,
      levels = c(
        "Direct",
        "Fix Direct",
        "FH",
        "Fixed FH",
        "GVF",
        "GVF-FH",
        "Joint"
      ),
      ordered = TRUE
    )
  )

admin1_medians <- unfixed_admin1$res.admin1 %>%
  dplyr::select(
    admin1.name,
    admin1_median = direct.est
  )

interval_width_comparison <- combined_data %>%
  dplyr::mutate(interval_width = upper - lower) %>%
  dplyr::group_by(model, area_type) %>%
  dplyr::summarise(
    mean_interval_width = mean(interval_width, na.rm = TRUE),
    median_interval_width = stats::median(interval_width, na.rm = TRUE),
    .groups = "drop"
  )

create_admin1_plot <- function(admin1_name) {
  plot_data <- combined_data %>%
    dplyr::filter(admin1.name == admin1_name) %>%
    dplyr::left_join(
      admin1_medians,
      by = "admin1.name"
    )
  
  current_median <- unique(plot_data$admin1_median)
  
  ggplot(
    plot_data,
    aes(x = reorder(admin2.name.full, as.numeric(area_type)))
  ) +
    geom_rect(
      data = plot_data %>%
        dplyr::distinct(area_type) %>%
        dplyr::mutate(
          xmin = -Inf,
          xmax = Inf,
          ymin = -Inf,
          ymax = Inf
        ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = area_type
      ),
      inherit.aes = FALSE,
      alpha = 0.15
    ) +
    geom_hline(
      yintercept = current_median,
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_pointrange(
      aes(
        y = estimate,
        ymin = lower,
        ymax = upper,
        color = model,
        shape = model
      ),
      position = position_dodge(width = 0.8),
      size = 0.5,
      fatten = 2
    ) +
    scale_color_manual(
      values = c(
        "Direct" = "#ff7f0e",
        "Fix Direct" = "#f768a1",
        "FH" = "#9467bd",
        "Fixed FH" = "#2ca02c",
        "GVF" = "#1f77b4",
        "GVF-FH" = "#8c564b",
        "Joint" = "#17becf"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Direct" = 17,
        "Fix Direct" = 18,
        "FH" = 3,
        "Fixed FH" = 15,
        "GVF" = 16,
        "GVF-FH" = 4,
        "Joint" = 8
      )
    ) +
    scale_fill_manual(
      values = c(
        "legal variance" = "white",
        "illegal variance" = "lightgrey",
        "missing" = "grey70"
      ),
      guide = "none"
    ) +
    labs(
      title = paste("Admin1:", admin1_name),
      subtitle = paste("Admin1 Median:", round(current_median, 4)),
      x = NULL,
      y = "Estimate with 95% CI",
      color = "Model",
      shape = "Model"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(color = "red", size = 10),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    facet_grid(
      . ~ area_type,
      scales = "free_x",
      space = "free_x"
    )
}

admin_list <- unique(combined_data$admin1.name)
admin1_interval_plots <- stats::setNames(
  lapply(admin_list, create_admin1_plot),
  admin_list
)

## **************************
## GVF plots ----------
## **************************
status_levels <- c(
  "Originally legal",
  "Originally illegal"
)

gvf_point_scatter_data <- gvf_point_scatter_data %>%
  dplyr::mutate(
    original_status = factor(
      original_status,
      levels = status_levels
    )
  )

gvf_variance_scatter_data <- gvf_variance_scatter_data %>%
  dplyr::mutate(
    original_status = factor(
      original_status,
      levels = status_levels
    )
  )

vsalm_point_scatter_data <- vsalm_point_scatter_data %>%
  dplyr::mutate(
    original_status = factor(
      original_status,
      levels = status_levels
    )
  )

vsalm_ci_width_scatter_data <- vsalm_ci_width_scatter_data %>%
  dplyr::mutate(
    original_status = factor(
      original_status,
      levels = status_levels
    )
  )

## Common x/y ranges
gvf_point_range <- range(
  c(
    gvf_point_scatter_data$fixed_estimate,
    gvf_point_scatter_data$gvf_estimate
  ),
  na.rm = TRUE
)

gvf_variance_range <- range(
  c(
    gvf_variance_scatter_data$fixed_logit_variance,
    gvf_variance_scatter_data$gvf_logit_variance
  ),
  na.rm = TRUE
)

vsalm_point_range <- range(
  c(
    vsalm_point_scatter_data$fixed_estimate,
    vsalm_point_scatter_data$vsalm_estimate
  ),
  na.rm = TRUE
)

vsalm_ci_width_range <- range(
  c(
    vsalm_ci_width_scatter_data$fixed_ci_width,
    vsalm_ci_width_scatter_data$vsalm_ci_width
  ),
  na.rm = TRUE
)

p_gvf_point <- ggplot(
  gvf_point_scatter_data,
  aes(
    x = fixed_estimate,
    y = gvf_estimate,
    color = original_status
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  coord_equal(
    xlim = gvf_point_range,
    ylim = gvf_point_range
  ) +
  labs(
    title = "GVF Point Estimate Comparison",
    subtitle = "Original probability scale",
    x = "Fixed direct estimate",
    y = "GVF estimate",
    color = "Original status"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

p_gvf_variance <- ggplot(
  gvf_variance_scatter_data,
  aes(
    x = fixed_logit_variance,
    y = gvf_logit_variance,
    color = original_status,
    size = n_i
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_point(
    alpha = 0.75
  ) +
  scale_size_continuous(
    range = c(1, 5)
  ) +
  coord_equal(
    xlim = gvf_variance_range,
    ylim = gvf_variance_range
  ) +
  labs(
    title = "GVF Logit-Scale Variance Comparison",
    subtitle = "Bubble size represents admin2 sample size",
    x = "Fixed direct logit variance",
    y = "GVF logit variance",
    color = "Original status",
    size = expression(n[i])
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

## **************************
## VSALM plots ----------
## **************************

p_vsalm_point <- ggplot(
  vsalm_point_scatter_data,
  aes(
    x = fixed_estimate,
    y = vsalm_estimate,
    color = original_status
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~method
  ) +
  coord_equal(
    xlim = vsalm_point_range,
    ylim = vsalm_point_range
  ) +
  labs(
    title = "VSALM Point Estimate Comparison",
    subtitle = "Original probability scale",
    x = "Fixed direct estimate",
    y = "VSALM posterior mean",
    color = "Original status"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

p_vsalm_ci_width <- ggplot(
  vsalm_ci_width_scatter_data,
  aes(
    x = fixed_ci_width,
    y = vsalm_ci_width,
    color = original_status
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~method
  ) +
  coord_equal(
    xlim = vsalm_ci_width_range,
    ylim = vsalm_ci_width_range
  ) +
  labs(
    title = "95% Interval Width Comparison",
    subtitle = paste0(
      "Fixed direct: 95% CI; ",
      "VSALM: 95% posterior credible interval"
    ),
    x = "Fixed direct 95% CI width",
    y = "VSALM 95% posterior credible interval width",
    color = "Original status"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

## **************************
## Print plots ----------
## **************************

print(p_gvf_point)
print(p_gvf_variance)

print(p_vsalm_point)
print(p_vsalm_ci_width)

for (admin_name in names(admin1_interval_plots)) {
  print(admin1_interval_plots[[admin_name]])
}

## Save Admin1 interval plots 
admin1_plot_dir <- here::here(
  "figures",
  "dhs",
  "new"
)

dir.create(
  admin1_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (admin_name in names(admin1_interval_plots)) {
  
  file_name <- paste0(
    "Nested_Admin1_",
    gsub("[^A-Za-z0-9]+", "_", admin_name),
    ".png"
  )
  
  ggplot2::ggsave(
    filename = file.path(
      admin1_plot_dir,
      file_name
    ),
    plot = admin1_interval_plots[[admin_name]],
    width = 12,
    height = 6,
    dpi = 300
  )
}

## **************************
## Save outputs ----------
## **************************

# output_dir <- here::here(
#   "figures",
#   "dhs",
#   "zambia_admin2_real_direct_gvf_vsalm"
# )
#
# dir.create(
#   output_dir,
#   recursive = TRUE,
#   showWarnings = FALSE
# )
#
# saveRDS(
#   list(
#     fixed_input = fixed_input,
#     direct_results = direct_results,
#
#     gvf_fit = gvf_fit,
#     gvf_results = gvf_results,
#     gvf_point_scatter_data = gvf_point_scatter_data,
#     gvf_variance_scatter_data = gvf_variance_scatter_data,
#
#     vsalm_bym2_fit = vsalm_bym2_fit,
#     vsalm_bym2_results = vsalm_bym2_results,
#     vsalm_results = vsalm_results,
#     vsalm_point_scatter_data = vsalm_point_scatter_data,
#     vsalm_ci_width_scatter_data = vsalm_ci_width_scatter_data,
#
#     method_results = method_results,
#     admin2_adj_mat = admin2_adj_mat,
#     diagnostics = analysis_diagnostics,
#
#     plots = list(
#       gvf_point = p_gvf_point,
#       gvf_variance = p_gvf_variance,
#       vsalm_point = p_vsalm_point,
#       vsalm_ci_width = p_vsalm_ci_width
#     )
#   ),
#   file.path(
#     output_dir,
#     "zambia_admin2_real_direct_gvf_vsalm_results.rds"
#   )
# )
#
# ggsave(
#   filename = file.path(
#     output_dir,
#     "gvf_point_comparison.png"
#   ),
#   plot = p_gvf_point,
#   width = 7,
#   height = 6,
#   dpi = 300
# )
#
# ggsave(
#   filename = file.path(
#     output_dir,
#     "gvf_logit_variance_comparison.png"
#   ),
#   plot = p_gvf_variance,
#   width = 8,
#   height = 6,
#   dpi = 300
# )
#
# ggsave(
#   filename = file.path(
#     output_dir,
#     "vsalm_point_comparison.png"
#   ),
#   plot = p_vsalm_point,
#   width = 10,
#   height = 6,
#   dpi = 300
# )
#
# ggsave(
#   filename = file.path(
#     output_dir,
#     "vsalm_ci_width_comparison.png"
#   ),
#   plot = p_vsalm_ci_width,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

