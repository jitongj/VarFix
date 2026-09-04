# *******************
# date: 2026.08.13
# task: This script is for direct-level GVF benchmark simulation of Zambia,
# author: Jitong Jiang
# ********************

library(dplyr)
library(purrr)
library(ggrepel)
library(devtools)
library(surveyPrev)
library(INLA)
library(geodata)
library(ggpattern)
library(SUMMER)
library(rdhs)
library(ggplot2)
library(patchwork)
library(tidyr)
library(kableExtra)
library(png)
library(grid) 
library(sf)
library(viridis)
library(gridExtra)
library(here)


## ======= sim_f1_s1_p1_r0.042_nationalWeight.RData ======= ## You can just load this data and run the plots!
## load(here::here("data/sim_f1_s1_p1_r0.042_nationalWeight.RData"))
scale_cluster_frame  <- 1.0    # eg 0.5, 1, 2
scale_cluster_sample <- 1.0    # eg 0.5, 1, 2
scale_population     <- 1.0    # eg 0.5, 1, 2
prevalence_value     <- 0.042  # prevalence rate
phantom_weight_quantile_value <- 0.10
## ================================= ##


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
poly.adm1 <- readRDS(here::here("data", "poly.adm1.rds"))
poly.adm2 <- readRDS(here::here("data", "poly.adm2.rds"))

# poly.adm1 <- geodata::gadm(country=country.abbrev, level=1, path=tempdir())
# poly.adm1 <- sf::st_as_sf(poly.adm1)
# poly.adm2 <- geodata::gadm(country=country.abbrev, level=2, path=tempdir())
# poly.adm2 <- sf::st_as_sf(poly.adm2) %>%
#   mutate(admin2.name.full = paste0(NAME_1, "_", NAME_2))
# # poly.adm2=poly.adm2[poly.adm2$ENGTYPE_2=="Local Authority",]



cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, 
                            by.adm1 = "NAME_1",by.adm2 = "NAME_2")
admin.info1 <- adminInfo(poly.adm = poly.adm1, admin = 1, by.adm = "NAME_1")
admin.info2 <- adminInfo(poly.adm = poly.adm2, admin = 2,
                         by.adm = "NAME_2", by.adm.upper = "NAME_1")


# **********************
# 2. Prepare pop den for admin1 and admin2   -------
# **********************
setwd(here::here("data", country))

country_shp_analysis <- readRDS('country_shp_analysis.rds')

pop.abbrev <- tolower(country.abbrev)
pop_file <- paste0(pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')

if(!file.exists(pop_file)){
  if (frame_year < 2021){
    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                  frame_year, "/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')
    
    download.file(url, pop_file, method = "libcurl", mode = "wb")
    
  } else{
    # download female and male tif
    url_female <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/",
                         frame_year, "/total_female_male/", toupper(pop.abbrev), "/",      
                         pop.abbrev,'_f_total_', frame_year,'_1km_UNadj.tif')
    
    url_male <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/",
                       frame_year, "/total_female_male/", toupper(pop.abbrev), "/",      
                       pop.abbrev,'_m_total_', frame_year,'_1km_UNadj.tif')
    
    file_female <- paste0("female_", pop.abbrev, "_", frame_year, ".tif")
    file_male <- paste0("male_", pop.abbrev, "_", frame_year, ".tif")
    
    download.file(url_female, file_female, method = "libcurl", mode = "wb")
    download.file(url_male, file_male, method = "libcurl", mode = "wb")
    
    # combine
    rast_female <- terra::rast(file_female)
    rast_male <- terra::rast(file_male)
    
    total_pop <- rast_female + rast_male
    terra::writeRaster(total_pop, pop_file, overwrite = TRUE)
  }
}

if (!exists("total_pop")) {
  total_pop <- terra::rast(pop_file)
}

#frame_pop <- terra::rast('nga_ppp_2006_1km_Aggregated_UNadj.tif') # pop at frame year
k0_1_pop <- total_pop #terra::rast('nga_k0_1_2018_1km.tif') # kid from 0 to 1 at survey year

pixel_grid <- as.data.frame(terra::xyFromCell(total_pop, 1:ncell(total_pop)))


get_pixel_adm_grid <- function(pixel_grid,
                               admin.level,
                               admin.poly){
  
  # Match stratification admin
  pixel_grid_sf <- st_as_sf(pixel_grid, coords = c("x", "y"),
                            crs = st_crs(admin.poly))
  
  adm_match <- st_join(pixel_grid_sf,admin.poly, join = st_intersects)
  
  pixel_grid_df <- data.frame(
    x = pixel_grid$x,
    y = pixel_grid$y
  )
  
  if(admin.level==1){
    pixel_grid_df$adm.name =adm_match[[paste0('NAME_1')]]
    pixel_grid_df$adm.name.full = pixel_grid_df$adm.name
  }else{
    
    pixel_grid_df$adm.name = adm_match[[paste0('NAME_',admin.level)]]
    pixel_grid_df$upper.adm.name = adm_match[[paste0('NAME_',admin.level-1)]]
    pixel_grid_df$adm.name.full = paste0(adm_match[[paste0('NAME_',admin.level-1)]],
                                         '_',
                                         adm_match[[paste0('NAME_',admin.level)]])
  }
  
  
  
  return(pixel_grid_df[complete.cases(pixel_grid_df),])
}

pixel_adm_grid <- get_pixel_adm_grid(pixel_grid = pixel_grid,
                                     admin.level = 2,
                                     admin.poly = country_shp_analysis[[2+1]])


# pixel pop
pixel_adm_grid$pop_den <- terra::extract(total_pop, pixel_adm_grid[,c('x','y')])[,2]
pixel_adm_grid$k0_1_pop <- terra::extract(k0_1_pop, pixel_adm_grid[,c('x','y')])[,2]
pixel_adm_grid <- pixel_adm_grid %>%
  filter(!is.na(pop_den))  


colnames(pixel_adm_grid) <- c("x","y","admin2.name","admin1.name", "admin2.name.full", "pop_den", "k0_1_pop")

# calculate admin1 total pop, k0-1
admin1_pop <- pixel_adm_grid %>%
  group_by(admin1.name) %>%
  summarise(
    total_pop_admin1 = sum(pop_den, na.rm = TRUE),
    k0_1_pop_admin1 = sum(k0_1_pop, na.rm = TRUE)
  )


# calculate admin2 total pop, k0-1
admin2_pop <- pixel_adm_grid %>%
  group_by(admin2.name.full) %>%
  summarise(
    total_pop_admin2 = sum(pop_den, na.rm = TRUE),
    k0_1_pop_admin2 = sum(k0_1_pop, na.rm = TRUE)
  )



# **********************
# 3. direct estimate   -------
# **********************


res_ad1 <- directEST_1030(data = data,
                     cluster.info = cluster.info,
                     admin = 1)
options(survey.adjust.domain.lonely = TRUE)
options(survey.lonely.psu = "adjust")
res_ad2 <- directEST_1030(data = data,
                     cluster.info = cluster.info,
                     admin = 2,
                     aggregation = FALSE,
                     var.fix = FALSE, all.fix=FALSE)
res_ad2_fix <- directEST_1030(data = data,
                                cluster.info = cluster.info,
                                admin = 2,
                                aggregation = FALSE,
                                var.fix = TRUE, all.fix=FALSE)

## individual level data
myData_tmp <- data0 %>%
  group_by(cluster) %>%   
  dplyr::summarise(      
    strata = unique(strata),  
    Ntrials = n(),      
    value = sum(value),   
    households_number = n_distinct(householdID)
  )
# combine individual level data with cluster info
myData <-merge(myData_tmp, cluster.info[["data"]], by="cluster")

# get admin2's direct estimate
direct <- res_ad2$res.admin2 %>%  
  filter(!is.na(direct.est)) %>% 
  distinct()


# **********************
# 4. simulation function   -------
# **********************

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

draw_cluster_sample_size <- function(
    mean_n = 17,
    sd_n = 5,
    min_n = 5,
    max_n = 30
) {
  n <- rnorm(1, mean_n, sd_n)
  while (n < min_n || n > max_n) {
    n <- rnorm(1, mean_n, sd_n)
  }
  round(n)
}

run_sampling_stratum_level <- function(
    admin1.name,
    cluster_frame,
    n_iterations = 1000,
    n_clusters = NULL,       # number of sampled cluster
    sigma0 = 0.5,            # Admin2 sigma
    sigma1 = 0.2,            # Cluster sigma
    sigma2 = 0.05,           # invidual sigma
    res_ad1 = NULL,          # Admin1 direct.est
    m = prevalence_value,
    coverage_level = 0.8,
    weight_target_n = 17     # planned/design sample size per cluster (must match
                             # draw_cluster_sample_size()'s mean_n); the per-person
                             # design weight is based on this fixed target, NOT on
                             # the realized (random) n_indiv_c -- otherwise the
                             # weight formula cancels out the realized count and
                             # every cluster ends up representing an (almost)
                             # identical population share, which is unrealistic.
) {
  

  cluster_sub <- cluster_frame[cluster_frame$admin1.name == admin1.name, ]
  cluster_ids <- cluster_sub$cluster_global_id
  M_c_total <- cluster_sub$M_c_total
  admin2_names_full <- cluster_sub$admin2.name.full
  admin2_names <- cluster_sub$admin2.name
  N <- nrow(cluster_sub)
  strata_vec <- cluster_sub$strata
  
 
  if (is.null(n_clusters)) {
    base_n <- sum(direct_merge$sampled_admin2_cluster[direct_merge$admin1.name == admin1.name], na.rm = TRUE)
   
    max_n <- admin2_cluster_pop %>%
      dplyr::filter(admin1.name == !!admin1.name) %>%
      dplyr::summarise(N = dplyr::n()) %>%
      dplyr::pull(N)
    
    n_clusters <- max(1L, min(as.integer(round(base_n * scale_cluster_sample)), max_n))
  }
  
  
  # baseline m
  if (is.null(m)) {
    m <- res_ad1$res.admin1$direct.est[res_ad1$res.admin1$admin1.name == admin1.name]
  }
  
  # alpha (admin2 level random effect) + e_c (cluster errir)
  alpha_admin2_map <- rnorm(length(unique(admin2_names_full)), 0, sigma0)
  names(alpha_admin2_map) <- unique(admin2_names_full)
  e_c <- rnorm(N, 0, sigma1)
  
 
  individual_df <- purrr::pmap_dfr(
    list(
      cluster_id = cluster_ids,
      M_c        = M_c_total,
      e_cluster  = e_c,
      admin2_full= admin2_names_full,
      admin2     = admin2_names,
      strata_c   = strata_vec          
    ),
    function(cluster_id, M_c, e_cluster, admin2_full, admin2, strata_c) { 
      alpha_admin2 <- alpha_admin2_map[admin2_full]
      e_ck <- rnorm(M_c, 0, sigma2)
      eta <- logit(m) + alpha_admin2 + e_cluster + e_ck   
      p_ck <- expit(eta)
      y <- rbinom(M_c, 1, p_ck)
      tibble(
        cluster_id = cluster_id,
        admin2.name.full = admin2_full,
        admin2.name = admin2,
        y = y,
        strata = strata_c            
      )
    }
  ) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::mutate(
      householdID = dplyr::row_number(),
      value = y,
      cluster = cluster_id,
      weight = 1,
      admin1.name = admin1.name
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(value, admin2.name.full, admin2.name, cluster, householdID, weight, admin1.name, strata)
  
  
  
  sampled_modt_list <- vector("list", n_iterations)
  
  for (i in 1:n_iterations) {
    # first stage PPS sampling cluster
    total_pop_admin1 <- unique(direct_merge$total_pop_admin1[direct_merge$admin1.name == admin1.name])
    sampled_clusters <- sample(
      cluster_ids,
      size = n_clusters,
      prob = n_clusters * M_c_total / total_pop_admin1,
      replace = FALSE
    )
    
    sampled_df <- individual_df[individual_df$cluster %in% sampled_clusters, ]
    
    # second stage: sample cluster-specific individuals and gives weight
    sampled_indiv <- sampled_df %>%
      group_by(cluster) %>%
      group_modify(~ {
        n_indiv_c <- draw_cluster_sample_size()
        .x %>%
          slice_sample(n = n_indiv_c) %>%
          mutate(n_indiv_c = n_indiv_c)
      }) %>%
      ungroup() %>%
      mutate(weight = round(total_pop_admin1 / (n_clusters * weight_target_n))) %>%
      select(-n_indiv_c)
    
    sampled_modt_list[[i]] <- sampled_indiv
  }
  
  return(list(
    full_population = individual_df,
    sampled_indiv_list = sampled_modt_list
  ))
}

# **********************
# 5. calculate some admin2 data   -------
# **********************

# Calculate the number of admin2 instances in each admin1 instance in dhs.
count_admin2 <- myData %>%
  group_by(admin1.name) %>%
  summarise(count_admin2 = length(unique(admin2.name.full)), .groups = "drop")

# Calculate the number of admin2 var << 1e-30 in each admin1 in dhs.
count_admin2_direct_var_0 <- direct %>%
  group_by(admin1.name) %>%
  summarise(count_admin2_direct_var_0 = sum(direct.var < 1e-30 , na.rm = TRUE), .groups = "drop")



# Calculate the proportion of admin2 with var = 0 in each admin1 in dhs.
count_admin2_direct_var_0 <- count_admin2_direct_var_0 %>%
  left_join(count_admin2, by = "admin1.name") %>%
  mutate(ratio = count_admin2_direct_var_0 / count_admin2)

# Calculate the number of clusters in each admin2 in dhs.
sampled_admin2_cluster <- myData %>%
  group_by(admin1.name, admin2.name, admin2.name.full) %>%
  summarise(sampled_admin2_cluster = n(), .groups = "drop")



# **********************
# 6. total number of clusters (i.e., in census)  -------
# **********************

#### Load admin1 total cluster:
frame_ea <- readRDS(here::here("data", country, paste0(pop.abbrev, "_frame_ea.rds")))


#### Estimate admin2 total cluster:

# Merge admin2's direct estimate
## Merge admin2's direct estimate
direct_merge <- merge(
  direct[, c("admin1.name", "admin2.name", "admin2.name.full",
             "direct.est", "direct.var","direct.se")], 
  count_admin2_direct_var_0, by = c("admin1.name")
)
direct_merge <- merge(
  direct_merge, sampled_admin2_cluster,
  by = c("admin1.name", "admin2.name", "admin2.name.full")
)
direct_merge <- merge(direct_merge, admin1_pop, by = "admin1.name")
direct_merge <- merge(direct_merge, admin2_pop, by = "admin2.name.full")

direct_merge <- direct_merge %>%
  dplyr::mutate(
    total_pop_admin1 = round(total_pop_admin1 * scale_population),
    total_pop_admin2 = round(total_pop_admin2 * scale_population)
  )


direct_merge <- direct_merge %>%
  dplyr::left_join(
    frame_ea %>% dplyr::select(admin1.name = strata, total_cluster_admin1 = total),
    by = "admin1.name"
  ) %>%
  dplyr::mutate(
    total_cluster_admin1_scaled = round(total_cluster_admin1 * scale_cluster_frame),
    
    ratio_admin2 = total_pop_admin2 / total_pop_admin1,
    estimated_cluster_admin2 = round(ratio_admin2 * total_cluster_admin1_scaled)
  )


poly.adm2 <- poly.adm2 %>%
  dplyr::mutate(admin2.name.full = paste(NAME_1, NAME_2, sep = "_"))


# **********************
# 7. number of individual per clusters (i.e., in census)   -------
# **********************
set.seed(2024)

admin2_cluster_pop <- purrr::pmap_dfr(
  list(
    admin1_name = direct_merge$admin1.name,
    admin2_name = direct_merge$admin2.name.full,
    total_pop = direct_merge$total_pop_admin2,
    N_clusters = direct_merge$estimated_cluster_admin2
  ),
  function(admin1_name, admin2_name, total_pop, N_clusters) {
    if (is.na(N_clusters) || N_clusters == 0) return(NULL)
    
    raw_weights <- rgamma(N_clusters, 1)
    cluster_weights <- raw_weights / sum(raw_weights)
    M_c_total <- pmax(round(cluster_weights * total_pop), 30)
    M_c_total <- round(M_c_total / sum(M_c_total) * total_pop)
    
    tibble(
      admin1.name = admin1_name,
      admin2.name.full = admin2_name,
      cluster_id_within_admin2 = 1:N_clusters,
      M_c_total = M_c_total
    )
  }
) %>%
  mutate(cluster_global_id = row_number())  



# **********************
# 8. direct merge   -------
# **********************

direct_merge_org <- direct_merge

direct_merge <- direct_merge %>%
  group_by(admin1.name) %>%
  mutate(
    mean_est_cluster = round(mean(estimated_cluster_admin2[estimated_cluster_admin2 >= 100], na.rm = TRUE)),
    estimated_cluster_admin2 = ifelse(estimated_cluster_admin2 < 100, mean_est_cluster, estimated_cluster_admin2),
    total_cluster_admin1 = sum(estimated_cluster_admin2, na.rm = TRUE),
    ratio_admin2 = estimated_cluster_admin2 / total_cluster_admin1,
    total_pop_admin1 = unique(total_pop_admin1),
    total_pop_admin2 = round(total_pop_admin1 * ratio_admin2),
    estimated_cluster_admin2 = round(estimated_cluster_admin2),
    total_cluster_admin1 = round(total_cluster_admin1)
  ) %>%
  ungroup() %>%
  select(-mean_est_cluster)





q_admin2 <- myData %>%
  dplyr::group_by(admin2.name.full) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_u = sum(strata == "urban", na.rm = TRUE),
    q_hat = ifelse(n > 0, n_u / n, NA_real_)
  ) %>%
  dplyr::ungroup()

admin2_cluster_pop <- purrr::pmap_dfr(
  list(
    admin1_name      = direct_merge$admin1.name,
    admin2_name_full = direct_merge$admin2.name.full,
    admin2_name      = direct_merge$admin2.name,   
    total_pop        = direct_merge$total_pop_admin2,
    N_clusters       = direct_merge$estimated_cluster_admin2
  ),
  function(admin1_name, admin2_name_full, admin2_name, total_pop, N_clusters) {
    if (is.na(N_clusters) || N_clusters == 0) return(NULL)
    
    # 1) Allocate cluster size
    raw_weights     <- rgamma(N_clusters, 1)
    cluster_weights <- raw_weights / sum(raw_weights)
    M_c_total <- pmax(round(cluster_weights * total_pop), 30L)
    M_c_total <- round(M_c_total / sum(M_c_total) * total_pop)
    
    # 2) City target percentage (if none, add it back to admin1 or the national logic).
    q2 <- q_admin2$q_hat[match(admin2_name_full, q_admin2$admin2.name.full)]
    q2 <- min(max(q2, 0), 1)
    
    # 3) Number of target city clusters
    N_u <- min(N_clusters, max(0, round(q2 * N_clusters)))
    
    # 4) Use a single stratum for simulated clusters.
    strata_vec <- rep("urban", N_clusters)
    
    tibble::tibble(
      admin1.name             = admin1_name,
      admin2.name.full        = admin2_name_full,
      admin2.name             = admin2_name,  
      cluster_id_within_admin2= seq_len(N_clusters),
      M_c_total               = as.integer(M_c_total),
      strata                  = strata_vec
    )
  }
) %>%
  dplyr::mutate(cluster_global_id = dplyr::row_number())


admin1_list <- unique(direct_merge$admin1.name)

# **********************
# 9. Run simulation for all admin2 -----
# **********************
process_simulation <- function(sampled_data, true_p_df) {
  # For this round of samples, count how many clusters this admin1 sampled from all admin2 and all strata.
  df_admin1 <- sampled_data %>%
    dplyr::group_by(admin1.name, strata) %>%
    dplyr::summarise(L_h = dplyr::n_distinct(cluster), .groups = "drop") %>%
    dplyr::group_by(admin1.name) %>%
    dplyr::summarise(
      K_admin1 = sum(L_h[L_h > 0], na.rm = TRUE),
      H_admin1 = sum(L_h > 0, na.rm = TRUE),
      df_unfixed_admin1 = as.integer(K_admin1 - H_admin1),
      df_fixed_admin1   = as.integer(K_admin1),
      .groups = "drop"
    )
  
  # 3 kinds direct results
  unfixed <- directEST_1030(
    data = sampled_data, NULL, admin = 2, aggregation = FALSE, var.fix = FALSE, all.fix = FALSE
  )$res.admin2
  
  fixed_triggered_obj <- directEST_1030(
    data = sampled_data, NULL, admin = 2,
    aggregation = FALSE,
    var.fix = TRUE,
    all.fix = FALSE,
    phantom_weight_quantile = phantom_weight_quantile_value
  )
  
  fixed_triggered <- fixed_triggered_obj$res.admin2
  fixed_area_names <- fixed_triggered_obj$fixed_areas
  
  fixed_all <- directEST_1030(
    data = sampled_data, NULL, admin = 2, aggregation = FALSE, var.fix = TRUE, all.fix = TRUE,
    phantom_weight_quantile = phantom_weight_quantile_value
  )$res.admin2

  area_sample_info <- sampled_data |>
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
  


  fixed_input <- fixed_triggered %>%
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
      )
    )

  # Add benchmark metadata without replacing any directEST-native result.
  benchmark_area_info <- unfixed_benchmark %>%
    dplyr::select(admin2.name.full, n_i, m_i, legal)

  proposed_all_fixed <- fixed_all %>%
    dplyr::left_join(benchmark_area_info, by = "admin2.name.full")

  proposed_illegal_fixed <- fixed_triggered %>%
    dplyr::left_join(benchmark_area_info, by = "admin2.name.full")

  # Sanity check:
  # all_fixed and fixed_illegal_unfixed_legal must keep the directEST-native
  # direct.est, direct.var, direct.logit.est, and direct.logit.var. They must
  # never be constructed as unfixed direct.est + fixed direct.var.
  
  keep_cols <- c(
    "admin2.name.full","admin1.name","admin2.name",
    "direct.est","direct.var","direct.var.gvf","direct.logit.est","direct.logit.var","n_i","m_i","legal",
    "prop_ph_denom","denom_phantom","denom_real","denom_all",
    "num_phantom","num_real","num_all","prop_ph_num",
    "p_hat_dataOnly","p_hat_withPh","delta_phantom"
  )
  
  res <- dplyr::bind_rows(
    unfixed_benchmark       %>% dplyr::mutate(type = "unfixed", method = "all_unfixed") %>% dplyr::select(dplyr::any_of(keep_cols), type, method),
    proposed_all_fixed      %>% dplyr::mutate(type = "fixed_all", method = "all_fixed") %>% dplyr::select(dplyr::any_of(keep_cols), type, method),
    proposed_illegal_fixed  %>% dplyr::mutate(type = "fixed_triggered", method = "fixed_illegal_unfixed_legal") %>% dplyr::select(dplyr::any_of(keep_cols), type, method),
    gvf_all                  %>% dplyr::mutate(type = "GVF", method = "GVF") %>% dplyr::select(dplyr::any_of(keep_cols), type, method)
  ) %>%
    dplyr::mutate(admin2.name.full = trimws(admin2.name.full)) %>%
    dplyr::left_join(df_admin1, by = "admin1.name") %>%
    dplyr::left_join(
      true_p_df %>% dplyr::transmute(admin2.name.full = trimws(admin2.name.full), true_p),
      by = "admin2.name.full"
    ) %>%
    dplyr::mutate(
      df_taylor = dplyr::case_when(
        type == "unfixed"         ~ df_unfixed_admin1,
        type == "fixed_triggered" ~ df_fixed_admin1,
        type == "fixed_all"       ~ df_fixed_admin1,
        type == "GVF"  ~ df_unfixed_admin1,
        TRUE ~ NA_integer_
      ),
      df_taylor = dplyr::if_else(is.finite(df_taylor) & df_taylor > 0, df_taylor, NA_integer_)
    )
  
  return(res)
}




run_simulations_and_collect_data <- function(admin1_list, admin2_cluster_pop, res_ad1) {
  set.seed(2024)
  
  sampled_by_admin1 <- purrr::map(admin1_list, function(ad1) {
    cat("Start admin1:", ad1, "...\n")
    
    res <- run_sampling_stratum_level(
      admin1.name   = ad1,
      cluster_frame = admin2_cluster_pop,
      n_iterations  = 1000,
      n_clusters    = NULL,
      res_ad1       = res_ad1,
      m             = prevalence_value
    )
    
    
    true_p_df <- res$full_population %>%
      dplyr::group_by(admin2.name.full, admin2.name, admin1.name) %>%
      dplyr::summarise(true_p = mean(value), .groups = "drop")
    cat("Finish admin1:", ad1, "\n")
    cat("---\n")
    
    list(
      res = res,
      true_p_df = true_p_df
    )
  })
  
  true_p_df <- sampled_by_admin1 %>%
    purrr::map("true_p_df") %>%
    dplyr::bind_rows()
  
  n_iterations <- length(sampled_by_admin1[[1]]$res$sampled_indiv_list)
  
  purrr::map_dfr(seq_len(n_iterations), function(sim_id) {
    if (sim_id %% 10 == 0) cat("Processing simulation:", sim_id, "/", n_iterations, "...\n")
    sampled_data <- sampled_by_admin1 %>%
      purrr::map_dfr(~ .x$res$sampled_indiv_list[[sim_id]])
    
    process_simulation(sampled_data, true_p_df) %>%
      dplyr::mutate(simulation_id = sim_id)
  })
}


simulation_start_time <- Sys.time()

all_detailed_results <- run_simulations_and_collect_data(
  admin1_list = admin1_list,
  admin2_cluster_pop = admin2_cluster_pop,
  res_ad1 = res_ad1
)

simulation_end_time <- Sys.time()
simulation_elapsed_seconds <- as.numeric(
  difftime(simulation_end_time, simulation_start_time, units = "secs")
)
cat(
  "Simulation runtime:\n",
  sprintf("  %.1f seconds\n", simulation_elapsed_seconds),
  sprintf("  %.2f minutes\n", simulation_elapsed_seconds / 60),
  sprintf("  %.2f hours\n", simulation_elapsed_seconds / 3600),
  sep = ""
)



calculate_performance_metrics <- function(all_results,
                                          conf_levels = c(0.5, 0.6, 0.70, 0.8, 0.95)) {
  
  base0 <- all_results %>%
    select(simulation_id, admin2.name.full, admin2.name, admin1.name,
           type, direct.est, direct.var, direct.logit.est, direct.logit.var,
           df_unfixed_admin1, df_fixed_admin1, legal, true_p)

  
  all_unfixed <- base0 %>%
    filter(type == "unfixed") %>%
    mutate(method = "all_unfixed", df_use = df_unfixed_admin1)

  all_fixed <- base0 %>%
    filter(type == "fixed_all") %>%
    mutate(method = "all_fixed", df_use = df_fixed_admin1)

  fixed_illegal_unfixed_legal <- base0 %>%
    filter(type == "fixed_triggered") %>%
    mutate(method = "fixed_illegal_unfixed_legal", df_use = df_fixed_admin1)

  gvf_benchmark <- base0 %>%
    filter(type == "GVF") %>%
    mutate(method = "GVF", df_use = df_unfixed_admin1)

  method_levels <- c(
    "all_unfixed", "all_fixed", "fixed_illegal_unfixed_legal",
    "GVF"
  )

  base_methods <- bind_rows(
    all_unfixed,
    all_fixed,
    fixed_illegal_unfixed_legal,
    gvf_benchmark
  ) %>%
    mutate(method = factor(method, levels = method_levels))
  
  
  one_level_summary <- function(df_in, conf_level) {
    alpha <- 1 - conf_level
    
    
    z_lo <- stats::qnorm((1 - conf_level) / 2)
    z_hi <- stats::qnorm(1 - (1 - conf_level) / 2)
    

    expit <- function(x) 1 / (1 + exp(-x))
    
    df <- df_in %>%
      mutate(
        # logit var and se
        var_logit_safe = tidyr::replace_na(direct.logit.var, 0),
        se_logit = sqrt(pmax(var_logit_safe, 0)),

        # Legacy point estimate for the original three methods.
        p_hat_point = expit(direct.logit.est),

        # A direct-estimate fallback is used only for an invalid GVF row.
        p_hat_point_gvf = if_else(
          is.finite(direct.logit.est),
          expit(direct.logit.est),
          direct.est
        ),
        
        # t distribution
        t_lo_raw = if_else(
          !is.na(df_use) & df_use > 0,
          stats::qt((1 - conf_level) / 2, df_use),
          z_lo
        ),
        t_hi_raw = if_else(
          !is.na(df_use) & df_use > 0,
          stats::qt(1 - (1 - conf_level) / 2, df_use),
          z_hi
        ),
        
        # Construct on the logit scale, then transform back to probability.
        ci_t_lb0 = expit(direct.logit.est + t_lo_raw * se_logit),
        ci_t_ub0 = expit(direct.logit.est + t_hi_raw * se_logit),
        ci_z_lb0 = expit(direct.logit.est + z_lo     * se_logit),
        ci_z_ub0 = expit(direct.logit.est + z_hi     * se_logit),
        
        # Match the legacy fallback for the first three methods: only an illegal
        # all_unfixed row is set to point estimate / zero-width interval, using
        # the existing legal column.  Failed/invalid Gao-style GVF rows receive
        # the same protection, without changing the other benchmark methods.
        gvf_invalid_ci = method == "GVF" &
          (!is.finite(direct.var) |
             direct.var < 0 |
             !is.finite(direct.logit.est) |
             !is.finite(direct.logit.var) |
             direct.logit.var < 0),
        use_point_interval = (method == "all_unfixed" &
                                (is.na(legal) | !legal)) |
          gvf_invalid_ci,
        point_interval_value = if_else(
          gvf_invalid_ci,
          p_hat_point_gvf,
          p_hat_point
        ),
        ci_t_lb = if_else(
          use_point_interval,
          point_interval_value,
          ci_t_lb0
        ),
        ci_t_ub = if_else(
          use_point_interval,
          point_interval_value,
          ci_t_ub0
        ),
        ci_z_lb = if_else(
          use_point_interval,
          point_interval_value,
          ci_z_lb0
        ),
        ci_z_ub = if_else(
          use_point_interval,
          point_interval_value,
          ci_z_ub0
        ),
        
        # truncated [0, 1]
        ci_t_lb = pmax(0, pmin(1, ci_t_lb)),
        ci_t_ub = pmax(0, pmin(1, ci_t_ub)),
        ci_z_lb = pmax(0, pmin(1, ci_z_lb)),
        ci_z_ub = pmax(0, pmin(1, ci_z_ub)),
        
        # coverage
        covered_t = (true_p >= ci_t_lb) & (true_p <= ci_t_ub),
        covered_z = (true_p >= ci_z_lb) & (true_p <= ci_z_ub),
        
        # width
        width_t = ci_t_ub - ci_t_lb,
        width_z = ci_z_ub - ci_z_lb,
        
        # interval score
        under_t = pmax(0, ci_t_lb - true_p),
        over_t  = pmax(0, true_p - ci_t_ub),
        score_t = width_t + (2 / alpha) * (under_t + over_t),
        
        under_z = pmax(0, ci_z_lb - true_p),
        over_z  = pmax(0, true_p - ci_z_ub),
        score_z = width_z + (2 / alpha) * (under_z + over_z),
        
        error = direct.est - true_p,
        sq_err = error^2
      )
    
    df %>%
      group_by(admin2.name.full, admin1.name, method) %>%
      summarise(
        true_p               = dplyr::first(true_p),
        avg_est              = mean(direct.est, na.rm = TRUE),
        # bias                 = mean(error, na.rm = TRUE),
        mse                  = mean(sq_err, na.rm = TRUE),
        rmse                 = sqrt(mse),
        coverage_rate_t      = mean(covered_t, na.rm = TRUE),
        coverage_rate_z      = mean(covered_z, na.rm = TRUE),
        avg_width_t          = mean(width_t, na.rm = TRUE),
        avg_width_z          = mean(width_z, na.rm = TRUE),
        avg_interval_score_t = mean(score_t, na.rm = TRUE),
        avg_interval_score_z = mean(score_z, na.rm = TRUE),
        avg_variance         = mean(direct.var, na.rm = TRUE),
        empirical_variance   = stats::var(direct.est, na.rm = TRUE),
        df_mean              = mean(df_use, na.rm = TRUE),
        n_sim                = dplyr::n(),
        .groups = "drop"
      ) %>%
      mutate(conf_level = conf_level)
  }
  
  
  
  out <- map_dfr(conf_levels, ~ one_level_summary(base_methods, .x))
  
  ## Legacy all_unfixed illegal rate: attach it only to all_unfixed.
  illegal_rate <- all_unfixed %>%
    group_by(admin2.name.full) %>%
    summarise(fail_pct_unfixed = mean(!legal, na.rm = TRUE) * 100,
              .groups = "drop") %>%
    mutate(method = factor("all_unfixed", levels = method_levels))

  ## Legacy benchmark base for variance diagnostics.
  emp_unfixed <- base_methods %>%
    filter(method == "all_unfixed") %>%
    group_by(admin2.name.full, admin1.name) %>%
    summarise(
      empirical_variance_unfixed = stats::var(direct.est, na.rm = TRUE),
      .groups = "drop"
    )
  
  out <- out %>%
    left_join(illegal_rate, by = c("admin2.name.full", "method")) %>%
    left_join(emp_unfixed, by = c("admin2.name.full", "admin1.name")) %>%
    mutate(
      # Match the legacy summary: fixed-method empirical variance and every
      # reported-variance comparison use the all_unfixed empirical variance.
      empirical_variance = if_else(
        method %in% c("all_fixed", "fixed_illegal_unfixed_legal"),
        empirical_variance_unfixed,
        empirical_variance
      ),
      variance_discrepancy = avg_variance - empirical_variance_unfixed
      # relative_variance_discrepancy = if_else(
      #   is.finite(empirical_variance_unfixed) & empirical_variance_unfixed > 0,
      #   variance_discrepancy / empirical_variance_unfixed,
      #   NA_real_)
    ) %>%
    arrange(admin2.name.full, method, conf_level) %>%
    select(
      admin2.name.full, admin1.name, method, conf_level,
      # Coverage rate (t / z)
      coverage_rate_t, coverage_rate_z,
      # CI width (t / z)
      avg_width_t, avg_width_z,
      # Interval score(t / z)
      avg_interval_score_t, avg_interval_score_z,
      # Error and variance
      mse, rmse, avg_est, true_p, avg_variance, empirical_variance,
      # bias,

      variance_discrepancy,
      # relative_variance_discrepancy,
      # t df
      df_mean,
      # Illegal rate (only all_unfixed has a value)
      fail_pct_unfixed,
      n_sim
    )
  
  return(out)
}




all_summary_metrics <- calculate_performance_metrics(
  all_results = all_detailed_results,
  conf_levels = c(0.50, 0.60, 0.70, 0.80, 0.95)
)

# **********************
# 10. plots   -------
# **********************

## plotting settings / labels / method levels
.method5_levels <- c("all_unfixed","all_fixed","fixed_illegal_unfixed_legal","GVF")
.method5_labels <- c(
  all_unfixed = "All—Unfixed",
  all_fixed   = "All—Fixed",
  fixed_illegal_unfixed_legal = "Fixed(Illegal)+Unfixed(Legal)",
  GVF = "GVF"
)
.method5_colors <- c(
  "all_unfixed" = "#ff7f0e",
  "all_fixed"   = "#1f77b4",
  "fixed_illegal_unfixed_legal" = "#2ca02c",
  "GVF" = "#9467bd"
)

scale_method5_fill <- function() {
  cols <- setNames(.method5_colors, .method5_labels)
  scale_fill_manual(
    values = cols,
    breaks = .method5_labels,
    labels = .method5_labels,
    name   = "Method",
    drop   = FALSE
  )
}


plot_admin1_conf_08 <- function(metrics,
                                y_col,
                                title,
                                ylab,
                                add_target_line = FALSE,
                                add_illegal_label = FALSE) {
  conf_target <- 0.80
  

  pd <- metrics %>%
    filter(conf_level == conf_target,
           method %in% .method5_levels) %>%
    mutate(
      method = factor(method, levels = .method5_levels, labels = .method5_labels)
    )
  

  admin1_levels <- sort(unique(pd$admin1.name))
  pd <- pd %>%
    mutate(admin1.name = factor(admin1.name, levels = admin1_levels))
  

  top_y <- pd %>%
    group_by(admin1.name) %>%
    summarise(
      y_top = max(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_top = ifelse(is.finite(y_top), y_top, 0),
      y_top = y_top + 0.02
    )
  
  if (add_illegal_label) {
    
    label_y <- 0.89
    
    annot <- pd %>%
      filter(method == "All—Unfixed") %>%
      group_by(admin1.name) %>%
      summarise(
        fail_rate = mean(fail_pct_unfixed, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        label = ifelse(
          is.na(fail_rate),
          "NA% illegal",
          sprintf("%.1f%% illegal", fail_rate)
        ),
        y_label = label_y
      )
    
  } else {
    annot <- NULL
  }
  

  p <- ggplot(
    pd,
    aes(x = admin1.name, y = .data[[y_col]], fill = method)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.7),
      width = 0.55,
      alpha = 0.95,
      outlier.shape = NA
    ) +
    scale_method5_fill() +
    labs(title = title, x = "Admin1 Region", y = ylab) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  

  if (add_illegal_label && !is.null(annot)) {
    p <- p + geom_text(
      data = annot,
      aes(
        x = admin1.name,
        y = y_label,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3,
      vjust = 1,
      color = "grey20"
    )
  }
  

  if (add_target_line) {
    p <- p + geom_hline(yintercept = conf_target,
                        linetype = "dashed",
                        color = "red")
  }
  
  p
}

## specialized plotting functions
plot_admin1_bias_var_est <- function(metrics) {
  p <- plot_admin1_conf_08(
    metrics,
    y_col = "variance_discrepancy",
    title = "Bias of Variance Estimators",
    ylab  = "Bias of Variance Estimators",
    add_target_line = FALSE,
    add_illegal_label = FALSE
  )
  p + geom_hline(yintercept = 0, linetype = "dashed")
}

plot_admin1_cov_normal_08 <- function(metrics) {
  p <- plot_admin1_conf_08(
    metrics,
    y_col = "coverage_rate_z",
    title = "Coverage",
    ylab  = "Empirical",
    add_target_line = TRUE,
    add_illegal_label = TRUE
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
    ylab  = "Average CI Width",
    add_illegal_label = FALSE
  )
}

plot_admin1_is_normal_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "avg_interval_score_z",
    title = "Interval Score",
    ylab  = "Interval Score",
    add_illegal_label = FALSE
  )
}

plot_admin1_cov_normal_095 <- function(metrics) {
  
  conf_target <- 0.95
  
  pd <- metrics %>%
    filter(
      conf_level == conf_target,
      method %in% .method5_levels
    ) %>%
    mutate(
      method = factor(
        method,
        levels = .method5_levels,
        labels = .method5_labels
      )
    )
  
  admin1_levels <- sort(unique(pd$admin1.name))
  
  pd <- pd %>%
    mutate(
      admin1.name = factor(
        admin1.name,
        levels = admin1_levels
      )
    )
  
  annot <- pd %>%
    filter(method == "All—Unfixed") %>%
    group_by(admin1.name) %>%
    summarise(
      fail_rate = mean(fail_pct_unfixed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = ifelse(
        is.na(fail_rate),
        "NA% illegal",
        sprintf("%.1f%% illegal", fail_rate)
      ),
      y_label = 1.00
    )
  
  p <- ggplot(
    pd,
    aes(
      x = admin1.name,
      y = coverage_rate_z,
      fill = method
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.7),
      width = 0.55,
      alpha = 0.95,
      outlier.shape = NA
    ) +
    geom_hline(
      yintercept = 0.95,
      linetype = "dashed",
      color = "red"
    ) +
    geom_text(
      data = annot,
      aes(
        x = admin1.name,
        y = y_label,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3,
      vjust = 1,
      color = "grey20"
    ) +
    scale_method5_fill() +
    labs(
      title = "Coverage",
      x = "Admin1 Region",
      y = "Empirical"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    coord_cartesian(ylim = c(0, 1.02))
  
  p
}

plot_admin1_width_normal_095 <- function(metrics) {
  
  conf_target <- 0.95
  
  pd <- metrics %>%
    filter(
      conf_level == conf_target,
      method %in% .method5_levels
    ) %>%
    mutate(
      method = factor(
        method,
        levels = .method5_levels,
        labels = .method5_labels
      )
    )
  
  admin1_levels <- sort(unique(pd$admin1.name))
  
  pd <- pd %>%
    mutate(
      admin1.name = factor(
        admin1.name,
        levels = admin1_levels
      )
    )
  
  ggplot(
    pd,
    aes(
      x = admin1.name,
      y = avg_width_z,
      fill = method
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.7),
      width = 0.55,
      alpha = 0.95,
      outlier.shape = NA
    ) +
    scale_method5_fill() +
    labs(
      title = "CI Width",
      x = "Admin1 Region",
      y = "Average CI Width"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

plot_admin1_is_normal_095 <- function(metrics) {
  
  conf_target <- 0.95
  
  pd <- metrics %>%
    filter(
      conf_level == conf_target,
      method %in% .method5_levels
    ) %>%
    mutate(
      method = factor(
        method,
        levels = .method5_levels,
        labels = .method5_labels
      )
    )
  
  admin1_levels <- sort(unique(pd$admin1.name))
  
  pd <- pd %>%
    mutate(
      admin1.name = factor(
        admin1.name,
        levels = admin1_levels
      )
    )
  
  ggplot(
    pd,
    aes(
      x = admin1.name,
      y = avg_interval_score_z,
      fill = method
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.7),
      width = 0.55,
      alpha = 0.95,
      outlier.shape = NA
    ) +
    scale_method5_fill() +
    labs(
      title = "Interval Score (95% CI)",
      x = "Admin1 Region",
      y = "Interval Score"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

## create plot objects
p_bias_var_est <- plot_admin1_bias_var_est(all_summary_metrics)
p_cov_z <- plot_admin1_cov_normal_08(all_summary_metrics)
p_wz    <- plot_admin1_width_normal_08(all_summary_metrics)
p_sz    <- plot_admin1_is_normal_08(all_summary_metrics)
p_cov_z_095 <- plot_admin1_cov_normal_095(all_summary_metrics)
p_wz_095 <- plot_admin1_width_normal_095(all_summary_metrics)
p_sz_095 <- plot_admin1_is_normal_095(all_summary_metrics)

## print plots
print(p_cov_z)
print(p_wz)
print(p_sz)
print(p_bias_var_est)
print(p_cov_z_095)
print(p_wz_095)
print(p_sz_095)


# **********************
# 11. save results   -------
# **********************

output_dir <- here::here("figures", "simulation_sensitivity", "nationalMean_10QWeight")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_output_dir <- here::here("data", "simulation_sensitivity")
dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  all_detailed_results,
  file.path(output_dir, "sensitivity10weight_simulation_direct_GVF_detailed_results.rds")
)
saveRDS(
  all_summary_metrics,
  file.path(output_dir, "sensitivity10weight_simulation_direct_GVF_performance_summary.rds")
)

# Keep data/simulation_sensitivity/sensitivity10weight_sim_f1_s1_p1_r0.042_nationalMean_10QWeight_GVF.RData (consumed by
# plot_stratified_legal_illegal_GVF_benchmark.R) in sync with this run.
save(
  all_detailed_results,
  all_summary_metrics,
  file = file.path(data_output_dir, "sensitivity10weight_sim_f1_s1_p1_r0.042_nationalMean_10QWeight_GVF.RData")
)

ggsave(
  filename = file.path(output_dir, "sensitivity10weight_simulation_direct_GVF_coverage.png"),
  plot = p_cov_z,
  width = 14, height = 6, dpi = 300
)
ggsave(
  filename = file.path(output_dir, "sensitivity10weight_simulation_direct_GVF_ci_width.png"),
  plot = p_wz,
  width = 14, height = 6, dpi = 300
)
ggsave(
  filename = file.path(output_dir, "sensitivity10weight_simulation_direct_GVF_interval_score.png"),
  plot = p_sz,
  width = 14, height = 6, dpi = 300
)
ggsave(
  filename = file.path(
    output_dir,
    "sensitivity10weight_simulation_direct_GVF_bias_variance_estimators.png"
  ),
  plot = p_bias_var_est,
  width = 14,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(
    output_dir,
    "sensitivity10weight_simulation_direct_GVF_coverage_95.png"
  ),
  plot = p_cov_z_095,
  width = 14,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(
    output_dir,
    "sensitivity10weight_simulation_direct_GVF_ci_width_95.png"
  ),
  plot = p_wz_095,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "sensitivity10weight_simulation_direct_GVF_interval_score_95.png"
  ),
  plot = p_sz_095,
  width = 14,
  height = 6,
  dpi = 300
)
