# *******************
# date: 2026.01.21
# task: This script is for simulation data of Zambia, 
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

run_sampling_stratum_level <- function(
    admin1.name,
    cluster_frame,
    n_iterations = 1000,
    n_clusters = NULL,       # number of sampled cluster
    n_indiv = 30,            # numbered of sampled individual in each cluster
    sigma0 = 0.5,            # Admin2 sigma
    sigma1 = 0.2,            # Cluster sigma
    sigma2 = 0.05,           # invidual sigma
    res_ad1 = NULL,          # Admin1 direct.est
    m = prevalence_value,
    coverage_level = 0.8,
    seed = 2024
) {
  set.seed(seed)
  

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
    
    # second stage: sample n_indiv in each cluster and gives weight
    sampled_indiv <- sampled_df %>%
      group_by(cluster) %>%
      slice_sample(n = n_indiv) %>%
      ungroup() %>%
      mutate(weight = round(total_pop_admin1 / (n_clusters * n_indiv)))
    
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
    
    # 4) Select city clusters by population weighting
    if (N_u == 0) {
      strata_vec <- rep("rural", N_clusters)
    } else if (N_u == N_clusters) {
      strata_vec <- rep("urban", N_clusters)
    } else {
      urban_idx  <- sample(seq_len(N_clusters), size = N_u, prob = M_c_total, replace = FALSE)
      strata_vec <- ifelse(seq_len(N_clusters) %in% urban_idx, "urban", "rural")
    }
    
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
  
  fixed_triggered <- directEST_1030(
    data = sampled_data, NULL, admin = 2, aggregation = FALSE, var.fix = TRUE, all.fix = FALSE
  )$res.admin2
  
  fixed_all <- directEST_1030(
    data = sampled_data, NULL, admin = 2, aggregation = FALSE, var.fix = TRUE, all.fix = TRUE
  )$res.admin2
  
  keep_cols <- c(
    "admin2.name.full","admin1.name","admin2.name",
    "direct.est","direct.var","direct.logit.est","direct.logit.var",
    "prop_ph_denom","denom_phantom","denom_real","denom_all",
    "num_phantom","num_real","num_all","prop_ph_num",
    "p_hat_dataOnly","p_hat_withPh","delta_phantom"
  )
  
  res <- dplyr::bind_rows(
    unfixed         %>% dplyr::mutate(type = "unfixed")          %>% dplyr::select(dplyr::any_of(keep_cols), type),
    fixed_triggered %>% dplyr::mutate(type = "fixed_triggered")  %>% dplyr::select(dplyr::any_of(keep_cols), type),
    fixed_all       %>% dplyr::mutate(type = "fixed_all")        %>% dplyr::select(dplyr::any_of(keep_cols), type)
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
        TRUE ~ NA_integer_
      ),
      df_taylor = dplyr::if_else(is.finite(df_taylor) & df_taylor > 0, df_taylor, NA_integer_)
    )
  
  return(res)
}




run_simulations_and_collect_data <- function(admin1_list, admin2_cluster_pop, res_ad1) {
  purrr::map_dfr(admin1_list, function(ad1) {
    cat("Start admin1:", ad1, "...\n")
    
    res <- run_sampling_stratum_level(
      admin1.name   = ad1,
      cluster_frame = admin2_cluster_pop,
      n_iterations  = 1000,
      n_clusters    = NULL,
      n_indiv       = 30,
      res_ad1       = res_ad1,
      m             = prevalence_value,
      seed          = 2024
    )
    
    
    true_p_df <- res$full_population %>%
      dplyr::group_by(admin2.name.full, admin2.name, admin1.name) %>%
      dplyr::summarise(true_p = mean(value), .groups = "drop")
    
    all_results <- purrr::map_dfr(
      res$sampled_indiv_list,
      ~process_simulation(.x, true_p_df),
      .id = "simulation"
    ) %>%
      dplyr::mutate(simulation_id = as.integer(simulation)) %>%
      dplyr::select(-simulation)
    
    cat("Finish admin1:", ad1, "\n")
    cat("---\n")
    all_results
  })
}

all_detailed_results <- run_simulations_and_collect_data(
  admin1_list = admin1_list,
  admin2_cluster_pop = admin2_cluster_pop,
  res_ad1 = res_ad1
)



calculate_performance_metrics <- function(all_results,
                                          conf_levels = c(0.5, 0.6, 0.70, 0.8, 0.90),
                                          legal_tol = 1e-12) {
  
  base0 <- all_results %>%
    mutate(legal = !is.na(direct.var) & (direct.var > legal_tol)) %>%
    select(simulation_id, admin2.name.full, admin2.name, admin1.name,
           type, direct.est, direct.var, direct.logit.est, direct.logit.var,
           df_unfixed_admin1, df_fixed_admin1, legal, true_p)
  
  ## 3 method types
  unfixed0         <- base0 %>% filter(type == "unfixed")
  fixed_all0       <- base0 %>% filter(type == "fixed_all")         
  fixed_triggered0 <- base0 %>% filter(type == "fixed_triggered")   
  
  if (nrow(fixed_all0) == 0) {
    stop("calculate_performance_metrics: No results found for type=='fixed_all'; please generate fixed_all upstream.")
  }
  if (nrow(fixed_triggered0) == 0) {
    stop("calculate_performance_metrics: No results found for type=='fixed_triggered'; please generate triggered results upstream.")
  }
  
  all_unfixed <- unfixed0 %>%
    mutate(method = "all_unfixed",
           df_use = df_unfixed_admin1)
  
  all_fixed <- fixed_all0 %>%
    mutate(method = "all_fixed",
           df_use = df_fixed_admin1)
  
  fixed_illegal_unfixed_legal <- fixed_triggered0 %>%
    mutate(method = "fixed_illegal_unfixed_legal",
           df_use = df_fixed_admin1)
  
  unfixed_legal <- unfixed0 %>%
    filter(legal) %>%
    mutate(method = "unfixed_legal",
           df_use = df_unfixed_admin1)
  
  base_methods <- bind_rows(
    all_unfixed %>% select(simulation_id, admin2.name.full, admin2.name, admin1.name,
                           method, direct.est, direct.var, direct.logit.est, direct.logit.var, df_use, legal, true_p),
    all_fixed   %>% select(simulation_id, admin2.name.full, admin2.name, admin1.name,
                           method, direct.est, direct.var, direct.logit.est, direct.logit.var, df_use, legal, true_p),
    fixed_illegal_unfixed_legal %>% select(simulation_id, admin2.name.full, admin2.name, admin1.name,
                                           method, direct.est, direct.var, direct.logit.est, direct.logit.var, df_use, legal, true_p),
    unfixed_legal %>% select(simulation_id, admin2.name.full, admin2.name, admin1.name,
                             method, direct.est, direct.var, direct.logit.est, direct.logit.var, df_use, legal, true_p)
  ) %>%
    mutate(method = factor(method,
                           levels = c("all_unfixed","all_fixed","fixed_illegal_unfixed_legal","unfixed_legal")))
  
  
  one_level_summary <- function(df_in, conf_level) {
    alpha <- 1 - conf_level
    
    
    z_lo <- stats::qnorm((1 - conf_level) / 2)
    z_hi <- stats::qnorm(1 - (1 - conf_level) / 2)
    

    expit <- function(x) 1 / (1 + exp(-x))
    
    df <- df_in %>%
      mutate(
        # logit var and se
        var_logit_safe = tidyr::replace_na(direct.logit.var, 0),
        se_logit       = sqrt(pmax(var_logit_safe, 0)),
        
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
        
        # The endpoints of the t/z interval at the logit scale and trans back
        ci_t_lb0 = expit(direct.logit.est + t_lo_raw * se_logit),
        ci_t_ub0 = expit(direct.logit.est + t_hi_raw * se_logit),
        ci_z_lb0 = expit(direct.logit.est + z_lo     * se_logit),
        ci_z_ub0 = expit(direct.logit.est + z_hi     * se_logit),
        
        
        p_hat_point = expit(direct.logit.est),
        
        # If it's all_unfixed and this is invalid, set the interval to point estimate to point estimate.
        ci_t_lb = if_else(
          method == "all_unfixed" & !legal,
          p_hat_point,
          ci_t_lb0
        ),
        ci_t_ub = if_else(
          method == "all_unfixed" & !legal,
          p_hat_point,
          ci_t_ub0
        ),
        ci_z_lb = if_else(
          method == "all_unfixed" & !legal,
          p_hat_point,
          ci_z_lb0
        ),
        ci_z_ub = if_else(
          method == "all_unfixed" & !legal,
          p_hat_point,
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
        
        # MSE
        sq_err = (direct.est - true_p)^2
      )
    
    df %>%
      group_by(admin2.name.full, admin1.name, method) %>%
      summarise(
        true_p               = dplyr::first(true_p),
        avg_est              = mean(direct.est, na.rm = TRUE),
        mse                  = mean(sq_err, na.rm = TRUE),
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
  
  ## all_unfixed illgeal rate
  illegal_rate <- all_unfixed %>%
    group_by(admin2.name.full) %>%
    summarise(fail_pct_unfixed = mean(!legal, na.rm = TRUE) * 100,
              .groups = "drop") %>%
    mutate(method = factor("all_unfixed",
                           levels = c("all_unfixed","all_fixed","fixed_illegal_unfixed_legal","unfixed_legal")))
  
  out <- out %>%
    left_join(illegal_rate, by = c("admin2.name.full","method"))
  
  ## Key: Forces the fixed empirical_variance to be equal to the all_unfixed empirical_variance
  emp_unfixed <- base_methods %>%
    filter(method == "all_unfixed") %>%
    group_by(admin2.name.full, admin1.name) %>%
    summarise(empirical_variance_unfixed = stats::var(direct.est, na.rm = TRUE),
              .groups = "drop")
  
  out <- out %>%
    left_join(emp_unfixed, by = c("admin2.name.full","admin1.name")) %>%
    mutate(
      empirical_variance = if_else(
        method %in% c("all_fixed","fixed_illegal_unfixed_legal"),
        empirical_variance_unfixed,
        empirical_variance
      ),
      
      variance_discrepancy = abs(empirical_variance_unfixed - avg_variance),
      relative_variance_discrepancy = if_else(
        !is.na(empirical_variance_unfixed) & empirical_variance_unfixed > 0,
        variance_discrepancy / empirical_variance_unfixed,
        NA_real_)
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
      mse, avg_est, true_p, avg_variance, empirical_variance,

      variance_discrepancy,
      relative_variance_discrepancy,
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
  conf_levels = c(0.50, 0.60, 0.70, 0.80, 0.90),
  legal_tol   = 1e-12
)


# all_detailed_results <- all_detailed_results %>%
#   dplyr::mutate(
#     admin1.name = dplyr::coalesce(.data[["admin1.name.x"]], .data[["admin1.name.y"]]),
#     admin2.name = dplyr::coalesce(.data[["admin2.name.x"]], .data[["admin2.name.y"]])
#   ) %>%
#   dplyr::select(-dplyr::matches("\\.x$|\\.y$"))
# 
# all_summary_metrics <- calculate_performance_metrics(
#   all_results = all_detailed_results,
#   conf_levels = c(0.70, 0.75, 0.80, 0.85, 0.90, 0.95),
#   legal_tol   = 1e-12
# )

# all_detailed_results %>%
#   dplyr::filter(type == "fixed") %>%
#   dplyr::select(simulation_id, admin1.name, admin2.name.full,
#                 K_raw, H, df_raw, phantom_used_urban, phantom_used_rural, df_taylor) %>%
#   head()
# 
# 
# df_summary = all_detailed_results %>%
#   dplyr::filter(type == "fixed") %>%
#   dplyr::select(simulation_id, admin1.name, admin2.name.full, df_taylor) 



# **********************
# 10. plots   -------
# **********************

.method3_levels <- c("all_unfixed","all_fixed","fixed_illegal_unfixed_legal")
.method3_labels <- c(
  all_unfixed = "All—Unfixed",
  all_fixed   = "All—Fixed",
  fixed_illegal_unfixed_legal = "Fixed(Illegal)+Unfixed(Legal)"
)
.method3_colors <- c(
  "all_unfixed" = "#ff7f0e",
  "all_fixed"   = "#1f77b4",
  "fixed_illegal_unfixed_legal" = "#2ca02c"
)

scale_method3_fill <- function() {
  cols <- setNames(.method3_colors, .method3_labels)
  scale_fill_manual(
    values = cols,
    breaks = .method3_labels,
    labels = .method3_labels,
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
           method %in% .method3_levels) %>%
    mutate(
      method = factor(method, levels = .method3_levels, labels = .method3_labels)
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
    annot <- pd %>%
      filter(method == "All—Unfixed") %>%
      group_by(admin1.name) %>%
      summarise(
        fail_rate = mean(fail_pct_unfixed, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        label = ifelse(is.na(fail_rate),
                       "NA% illegal",
                       sprintf("%.1f%% illegal", fail_rate))
      ) %>%
      left_join(top_y, by = "admin1.name")
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
    scale_method3_fill() +
    labs(title = title, x = "Admin1 Region", y = ylab) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  

  if (add_illegal_label && !is.null(annot)) {
    p <- p + geom_text(
      data = annot,
      aes(x = admin1.name, y = y_top, label = label),
      inherit.aes = FALSE,
      size = 3,
      vjust = 0,
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




# 1) Variance discrepancy
plot_admin1_var_disc_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "variance_discrepancy",
    title = "Variance Discrepancy (conf = 0.80)",
    ylab  = "Empirical − Reported Variance",
    add_target_line = FALSE,
    add_illegal_label = FALSE
  )
}

# 2) Coverage (Normal)
plot_admin1_cov_normal_08 <- function(metrics) {
  p <- plot_admin1_conf_08(
    metrics,
    y_col = "coverage_rate_z",
    title = "Coverage",
    ylab  = "Empirical",
    add_target_line = TRUE,
    add_illegal_label = TRUE
  )
  # p + coord_cartesian(ylim = c(0, 1))
  yvals <- metrics$coverage_rate_z
  y_min <- min(yvals, na.rm = TRUE)
  y_max <- max(yvals, na.rm = TRUE)
  
  p + coord_cartesian(ylim = c(y_min, y_max))
}

# 3) Coverage (t)
plot_admin1_cov_t_08 <- function(metrics) {
  p <- plot_admin1_conf_08(
    metrics,
    y_col = "coverage_rate_t",
    title = "Coverage (t) (conf = 0.80)",
    ylab  = "Empirical Coverage (t)",
    add_target_line = TRUE,
    add_illegal_label = TRUE
  )
  p + coord_cartesian(ylim = c(0, 1))
}

# 4) df (t)
plot_admin1_df_t_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "df_mean",
    title = "Degrees of Freedom (t) (conf = 0.80)",
    ylab  = "Mean df (per domain)",
    add_illegal_label = FALSE
  )
}

# 5) CI width (Normal)
plot_admin1_width_normal_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "avg_width_z",
    title = "CI Width",
    ylab  = "Average CI Width",
    add_illegal_label = FALSE
  )
}

# 6) CI width (t)
plot_admin1_width_t_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "avg_width_t",
    title = "CI Width (t) (conf = 0.80)",
    ylab  = "Average CI Width (t)",
    add_illegal_label = FALSE
  )
}

# 7) Interval score (Normal)
plot_admin1_is_normal_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "avg_interval_score_z",
    title = "Interval Score",
    ylab  = "Interval Score",
    add_illegal_label = FALSE
  )
}

# 8) Interval score (t)
plot_admin1_is_t_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "avg_interval_score_t",
    title = "Interval Score (t) (conf = 0.80)",
    ylab  = "Interval Score (t)",
    add_illegal_label = FALSE
  )
}

plot_admin1_relative_var_disc_08 <- function(metrics) {
  plot_admin1_conf_08(
    metrics,
    y_col = "relative_variance_discrepancy",
    title = "Relative Variance Discrepancy (conf = 0.80)",
    ylab  = "|Empirical − Reported| / Empirical Variance",
    add_target_line = FALSE,
    add_illegal_label = FALSE
  )
}



p_var   <- plot_admin1_var_disc_08(all_summary_metrics)
p_cov_z <- plot_admin1_cov_normal_08(all_summary_metrics)
p_cov_t <- plot_admin1_cov_t_08(all_summary_metrics)
p_df    <- plot_admin1_df_t_08(all_summary_metrics)
p_wz    <- plot_admin1_width_normal_08(all_summary_metrics)
p_wt    <- plot_admin1_width_t_08(all_summary_metrics)
p_sz    <- plot_admin1_is_normal_08(all_summary_metrics)
p_st    <- plot_admin1_is_t_08(all_summary_metrics)
p_rel_var <- plot_admin1_relative_var_disc_08(all_summary_metrics)


# print(p_var); print(p_cov_z); print(p_cov_t); print(p_df)
# print(p_wz);  print(p_wt);  print(p_sz);    print(p_st)
# print(p_rel_var)


### Here this the final plots!!!!!
print(p_cov_z)
print(p_wz)
print(p_sz)

## these mimic plots are saved under figures/simulation


