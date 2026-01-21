# *******************
# date: 2026.01.21
# task: This script is for example of Lavushimanda, national aggreagted estimate, figure7-10, and table2
# author: Jitong Jiang
# ********************

library(surveyPrev)
library(SUMMER)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sf)
library(tidyr)
library(stringr)
library(ggrepel)
library(tibble)


source(here::here("codes", "directEST_1030_national.R"))
source(here::here("codes", "fhModel_1030.R"))


# **********************
# 1. Basic setting   -------
# **********************
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


cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, 
                            by.adm1 = "NAME_1",by.adm2 = "NAME_2")
admin.info1 <- adminInfo(poly.adm = poly.adm1, admin = 1, by.adm = "NAME_1")
admin.info2 <- adminInfo(poly.adm = poly.adm2, admin = 2,
                         by.adm = "NAME_2", by.adm.upper = "NAME_1")

# **********************
# 2. Direct estimate results-------
# **********************
res_ad1 <- directEST_1030(data = data,
                          cluster.info = cluster.info,
                          admin = 1,
                          aggregation = FALSE,
                          var.fix = FALSE)

res_ad2 <- directEST_1030(data = data,
                          cluster.info = cluster.info,
                          admin = 2,
                          aggregation = FALSE,
                          var.fix = FALSE)

res_ad2_fix <- directEST_1030(data = data,
                              cluster.info = cluster.info,
                              admin = 2,
                              aggregation = FALSE,
                              var.fix = TRUE)

data_urban <- data[data$strata == "urban", ]
data_rural <- data[data$strata == "rural", ]

res_ad0 <- directEST(data = data,
                     cluster.info = cluster.info,
                     admin = 0)
res_ad0_ubran <- directEST(data = data_urban,
                           cluster.info = cluster.info,
                           admin = 0)
res_ad0_rural <- directEST(data = data_rural,
                           cluster.info = cluster.info,
                           admin = 0)


# **********************
# 3. FH estimate results-------
# **********************
bad_admin2 <- res_ad2_fix$fixed_areas
bad_clusters <- subset(cluster.info$data, admin2.name.full %in% bad_admin2)$cluster

smth_res_ad2 <- fhModel_1030(subset(data, !cluster %in% bad_clusters),
                             cluster.info = cluster.info,
                             admin.info = admin.info2,
                             admin = 2,
                             model = "bym2",
                             aggregation = FALSE,
                             var.fix = FALSE,
                             nested = TRUE)

smth_res_ad2_NoNested <- fhModel_1030(subset(data, !cluster %in% bad_clusters),
                                      cluster.info = cluster.info,
                                      admin.info = admin.info2,
                                      admin = 2,
                                      model = "bym2",
                                      aggregation = FALSE,
                                      var.fix = FALSE,
                                      nested = FALSE)


smth_res_ad2_fix <- fhModel_1030(data,
                                 cluster.info = cluster.info,
                                 admin.info = admin.info2,
                                 admin = 2,
                                 model = "bym2",
                                 aggregation = FALSE,
                                 var.fix = TRUE,
                                 nested = TRUE)

smth_res_ad2_fix_NoNested <- fhModel_1030(data,
                                          cluster.info = cluster.info,
                                          admin.info = admin.info2,
                                          admin = 2,
                                          model = "bym2",
                                          aggregation = FALSE,
                                          var.fix = TRUE,
                                          nested = FALSE)

# **********************
# 4. Classification by type -------
# **********************
country_shp_analysis <- readRDS(here::here("data", "Zambia", "country_shp_analysis.rds"))

data0 <- data %>% filter(!is.na(value))

myData_tmp <- data0 %>%
  group_by(cluster) %>%   # Group by cluster
  dplyr::summarise(       # Compute within-cluster summary statistics
    strata = unique(strata),  # Unique strata value per cluster (assuming one strata per cluster)
    Ntrials = n(),            # Number of individuals in the cluster
    value = sum(value),       # Total number of successes (sum of individual values)
    households_number = n_distinct(householdID) # Number of distinct households in the cluster
  )

# combine individual level data with cluster info
myData <-merge(myData_tmp, cluster.info[["data"]], by="cluster")

sampled_admin2_cluster <- myData %>%
  group_by(admin1.name, admin2.name, admin2.name.full) %>%
  summarise(sampled_admin2_cluster = n(), .groups = "drop")


full_admin2 <- country_shp_analysis[["Admin-2"]] |>
  dplyr::mutate(admin2.name.full = paste(NAME_1, NAME_2, sep = "_"))

unsampled_admin2 <- dplyr::filter(full_admin2, !(admin2.name.full %in% sampled_admin2_cluster$admin2.name.full))

sparse_admin2 <- bad_admin2

missing_admin2 <- unsampled_admin2$admin2.name.full



# **********************
# 5. Lavushimanda example   -------
# **********************

lavushimanda_ad2 <- bind_rows(
  smth_res_ad2$res.admin2 %>%
    filter(admin2.name.full == "Muchinga_Lavushimanda") %>%
    mutate(var.fix = FALSE, nested = TRUE),
  
  smth_res_ad2_NoNested$res.admin2 %>%
    filter(admin2.name.full == "Muchinga_Lavushimanda") %>%
    mutate(var.fix = FALSE, nested = FALSE),
  
  smth_res_ad2_fix$res.admin2 %>%
    filter(admin2.name.full == "Muchinga_Lavushimanda") %>%
    mutate(var.fix = TRUE, nested = TRUE),
  
  smth_res_ad2_fix_NoNested$res.admin2 %>%
    filter(admin2.name.full == "Muchinga_Lavushimanda") %>%
    mutate(var.fix = TRUE, nested = FALSE)
) %>%
  select(var.fix, nested, everything())



# **********************
# 6. Figure7 (Nested vs Non-Nested) -------
# **********************

## 6.1 Build comparison data frames for figure7 ----


# Non-fixed: subset data, nested vs non-nested
df_nonfix <- smth_res_ad2$res.admin2 %>%
  select(admin1.name, admin2.name.full, mean_nested = mean)

df_nonfix_nn <- smth_res_ad2_NoNested$res.admin2 %>%
  select(admin1.name, admin2.name.full, mean_nonnested = mean)

df_nonfix <- df_nonfix %>%
  inner_join(df_nonfix_nn,
             by = c("admin1.name", "admin2.name.full"))

# Fixed: full data, nested vs non-nested
df_fix <- smth_res_ad2_fix$res.admin2 %>%
  select(admin1.name, admin2.name.full, mean_nested = mean)

df_fix_nn <- smth_res_ad2_fix_NoNested$res.admin2 %>%
  select(admin1.name, admin2.name.full, mean_nonnested = mean)

df_fix <- df_fix %>%
  inner_join(df_fix_nn,
             by = c("admin1.name", "admin2.name.full"))



## 6.2. Plot figure7  -------


all_lims <- range(
  c(df_nonfix$mean_nested,
    df_nonfix$mean_nonnested,
    df_fix$mean_nested,
    df_fix$mean_nonnested),
  na.rm = TRUE
)


base_theme <- theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none"   # turn OFF legend here
  )

p_nonfix <- ggplot(df_nonfix,
                   aes(x = mean_nonnested,
                       y = mean_nested,
                       color = admin1.name)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = all_lims, ylim = all_lims, expand = FALSE) +
  labs(
    x = "Non-nested mean",
    y = "Nested mean",
    title = "Non-fixed variance"
  ) +
  base_theme

p_fix <- ggplot(df_fix,
                aes(x = mean_nonnested,
                    y = mean_nested,
                    color = admin1.name)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = all_lims, ylim = all_lims, expand = FALSE) +
  labs(
    x = "Non-nested mean",
    y = "Nested mean",
    title = "Fixed variance"
  ) +
  base_theme


nested_vs_nonnested <-
  (p_nonfix | p_fix) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Admin-2: Nested vs Non-nested Fay-Herriot Model") &
  theme(legend.position = "right")


ggsave(
  filename = here::here("figures", "dhs", "nested_vs_nonnested.png"),
  plot = nested_vs_nonnested,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)


# **********************
# 7. Figure8  -------
# **********************

## 7.1 Built dataframe for figure8  --------

full_admin2 <- country_shp_analysis[["Admin-2"]] |>
  dplyr::mutate(admin2.name.full = paste(NAME_1, NAME_2, sep = "_"))


full_admin2_nogeo <- full_admin2 %>%
  sf::st_drop_geometry() %>%
  select(admin2.name.full)

point_type_df <- full_admin2_nogeo %>%
  mutate(
    point_type = dplyr::case_when(
      admin2.name.full %in% missing_admin2 ~ "no data",
      admin2.name.full %in% sparse_admin2  ~ "illegal variance",
      TRUE                                 ~ "legal variance"
    )
  ) %>%
  mutate(
    point_type = factor(
      point_type,
      levels = c("legal variance", "no data", "illegal variance"),
      labels = c("legal variance", "No data", "Illegal variance")
    )
  )


## 7.2 FH mean vs Direct mean  --------

fh_direct_df <- smth_res_ad2$res.admin2 %>%
  select(admin2.name.full,
         fh_mean = mean) %>%
  full_join(  
    res_ad2$res.admin2 %>%
      select(admin2.name.full,
             direct_mean = direct.est),
    by = "admin2.name.full"
  ) %>%
  left_join(point_type_df, by = "admin2.name.full")


lims1 <- range(c(fh_direct_df$direct_mean,
                 fh_direct_df$fh_mean),
               na.rm = TRUE)

## 7.3 Fixed FH mean vs FH mean ---------
fh_fix_mean_df <- smth_res_ad2$res.admin2 %>%
  select(admin2.name.full,
         fh_mean = mean) %>%
  inner_join(
    smth_res_ad2_fix$res.admin2 %>%
      select(admin2.name.full,
             fh_fix_mean = mean),
    by = "admin2.name.full"
  ) %>%
  left_join(point_type_df, by = "admin2.name.full")

lims2 <- range(c(fh_fix_mean_df$fh_mean,
                 fh_fix_mean_df$fh_fix_mean),
               na.rm = TRUE)

## 7.4 Fixed FH SD vs FH SD ---------

fh_fix_sd_df <- smth_res_ad2$res.admin2 %>%
  select(admin2.name.full,
         fh_sd = sd) %>%
  inner_join(
    smth_res_ad2_fix$res.admin2 %>%
      select(admin2.name.full,
             fh_fix_sd = sd),
    by = "admin2.name.full"
  ) %>%
  left_join(point_type_df, by = "admin2.name.full")

lims3 <- range(c(fh_fix_sd_df$fh_sd,
                 fh_fix_sd_df$fh_fix_sd),
               na.rm = TRUE)

## 7.5 Plot figure8 ---------


base_theme <- theme_bw() +
  theme(
    aspect.ratio   = 1,
    legend.position = "none"
  )

# FH vs Direct
p1 <- ggplot(fh_direct_df,
             aes(x = direct_mean,
                 y = fh_mean,
                 color = point_type,
                 shape = point_type)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = lims1, ylim = lims1, expand = FALSE) +
  labs(
    x = "Direct estimate",
    y = "FH estimate",
    title = "FH vs Direct (means)"
  ) +
  base_theme

# Fixed FH mean vs FH mean
p2 <- ggplot(fh_fix_mean_df,
             aes(x = fh_mean,
                 y = fh_fix_mean,
                 color = point_type,
                 shape = point_type)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = lims2, ylim = lims2, expand = FALSE) +
  labs(
    x = "FH mean",
    y = "Fixed FH mean",
    title = "Fixed FH vs FH (means)"
  ) +
  base_theme

# Fixed FH SD vs FH SD
p3 <- ggplot(fh_fix_sd_df,
             aes(x = fh_sd,
                 y = fh_fix_sd,
                 color = point_type,
                 shape = point_type)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = lims3, ylim = lims3, expand = FALSE) +
  labs(
    x = "FH posterior SD",
    y = "Fixed FH posterior SD",
    title = "Fixed FH vs FH (SDs)"
  ) +
  base_theme


scatter_compare <-
  (p1 | p2 | p3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


ggsave(
  filename = here::here("figures", "dhs", "scatter_compare.png"),
  plot = scatter_compare,
  width = 12,
  height = 4,
  dpi = 300,
  bg = "white"
)



# **********************
# 8. Figure9  -------
# **********************

## 8.1 Raster pop data and aggregate--------

k0_5_pop <- terra::rast(
  here::here("data", "subpop", "zmb_k0_5_2018_1km.tif")
)

pixel_grid <- as.data.frame(terra::xyFromCell(k0_5_pop, 1:terra::ncell(k0_5_pop)))

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
pixel_adm_grid$pop_den <- terra::extract(k0_5_pop, pixel_adm_grid[,c('x','y')])[,2]
pixel_adm_grid <- pixel_adm_grid %>%
  filter(!is.na(pop_den)) 

colnames(pixel_adm_grid) <- c("x","y","admin2.name","admin1.name", "admin2.name.full", "pop_den")


## aggregate pixel pop

# calculate admin1 k0-5
admin1_pop <- pixel_adm_grid %>%
  group_by(admin1.name) %>%
  summarise(
    total_pop_admin1 = sum(pop_den, na.rm = TRUE) # total_pop is k0_5
  )

# calculate admin2 k0-5
admin2_pop <- pixel_adm_grid %>%
  group_by(admin2.name.full, admin1.name) %>%
  summarise(
    total_pop_admin2 = sum(pop_den, na.rm = TRUE) # total_pop is k0_5
  )


## 8.2 Aggregate model' result--------

aggregate_admin2_to_admin1 <- function(smth_res_ad2, admin2_pop,
                                       post_name = "admin2_post",
                                       res_name  = "res.admin2",
                                       admin2_key = "admin2.name.full",
                                       admin1_key = "admin1.name",
                                       pop_key    = "total_pop_admin2",
                                       probs_ci   = c(0.025, 0.975),
                                       tol_wsum   = 1e-10,
                                       set_colnames_if_missing = TRUE,
                                       reorder_post_if_needed  = TRUE,
                                       check_dims = TRUE,
                                       verbose = TRUE) {
  # Dependencies
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  `%>%` <- dplyr::`%>%`
  

  # 1) Extract posterior draws and admin2 summary (for target order)

  if (!post_name %in% names(smth_res_ad2)) stop("Missing smth_res_ad2[['", post_name, "']].")
  if (!res_name  %in% names(smth_res_ad2)) stop("Missing smth_res_ad2[['", res_name,  "']].")
  
  admin2_post <- smth_res_ad2[[post_name]]
  res_ad2     <- smth_res_ad2[[res_name]]
  
  if (!is.matrix(admin2_post) && !is.data.frame(admin2_post)) {
    stop("smth_res_ad2[['", post_name, "']] must be a matrix or data.frame.")
  }
  admin2_post <- as.matrix(admin2_post)
  
  if (!is.data.frame(res_ad2)) res_ad2 <- as.data.frame(res_ad2)
  
  if (!admin2_key %in% names(res_ad2)) {
    stop("Column '", admin2_key, "' not found in smth_res_ad2[['", res_name, "']].")
  }
  target_order <- res_ad2[[admin2_key]]
  
  if (check_dims) {
    if (ncol(admin2_post) != length(target_order)) {
      stop("Dimension mismatch: ncol(admin2_post) = ", ncol(admin2_post),
           " but length(target_order) = ", length(target_order), ".")
    }
  }
  
  # 2) Ensure admin2_post colnames match target order (or set them)

  if (is.null(colnames(admin2_post)) && set_colnames_if_missing) {
    colnames(admin2_post) <- target_order
  }
  
  if (reorder_post_if_needed && !identical(colnames(admin2_post), target_order)) {
    if (verbose) {
      warning("admin2_post colnames are not identical to res.admin2 order; reordering admin2_post to match res.admin2.")
    }
    idx_post <- match(target_order, colnames(admin2_post))
    if (anyNA(idx_post)) {
      missing_names <- target_order[is.na(idx_post)]
      stop("These admin2 are missing in admin2_post colnames:\n", paste(missing_names, collapse = ", "))
    }
    admin2_post <- admin2_post[, idx_post, drop = FALSE]
  }
  
  # 3) Reorder admin2_pop to match target order
  needed_cols <- c(admin2_key, admin1_key, pop_key)
  missing_cols <- setdiff(needed_cols, names(admin2_pop))
  if (length(missing_cols) > 0) {
    stop("admin2_pop is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  idx_pop <- match(target_order, admin2_pop[[admin2_key]])
  if (anyNA(idx_pop)) {
    missing_pop <- target_order[is.na(idx_pop)]
    stop("These admin2 are missing in admin2_pop:\n", paste(missing_pop, collapse = ", "))
  }
  
  admin2_pop_ord <- admin2_pop %>%
    dplyr::slice(idx_pop) %>%
    dplyr::mutate(!!admin2_key := factor(.data[[admin2_key]], levels = target_order))
  
  # Extra safety check
  if (!all(as.character(admin2_pop_ord[[admin2_key]]) == target_order)) {
    stop("Reordering failed: admin2_pop_ord is not perfectly aligned with target_order.")
  }
  
  
  # 4) Compute within-admin1 population weights (sum to 1 within each admin1)
  admin2_w <- admin2_pop_ord %>%
    dplyr::group_by(.data[[admin1_key]]) %>%
    dplyr::mutate(
      pop_admin1 = sum(.data[[pop_key]], na.rm = TRUE),
      weight     = .data[[pop_key]] / pop_admin1
    ) %>%
    dplyr::ungroup()
  
  check_w <- admin2_w %>%
    dplyr::group_by(.data[[admin1_key]]) %>%
    dplyr::summarise(w_sum = sum(weight), .groups = "drop")
  
  if (any(abs(check_w$w_sum - 1) > tol_wsum)) {
    if (verbose) print(check_w)
    stop("Within-admin1 weights do not sum to 1 (tolerance = ", tol_wsum, "). Check population inputs.")
  }
  
  # Preserve admin1 order by first appearance in the aligned admin2 list
  admin1_names <- admin2_w %>%
    dplyr::distinct(.data[[admin1_key]]) %>%
    dplyr::pull(.data[[admin1_key]])
  

  # 5) Aggregate admin2 posterior draws to admin1 using weights
  admin1_post <- sapply(admin1_names, function(a1) {
    idx <- which(admin2_w[[admin1_key]] == a1)
    w   <- admin2_w$weight[idx]
    as.numeric(admin2_post[, idx, drop = FALSE] %*% w)
  })
  
  admin1_post <- as.matrix(admin1_post)
  colnames(admin1_post) <- admin1_names
  rownames(admin1_post) <- rownames(admin2_post)
  
  # 6) Admin1 posterior summary
  p_lo <- probs_ci[1]
  p_hi <- probs_ci[2]
  
  admin1_summary <- data.frame(
    admin1.name = admin1_names,
    mean   = colMeans(admin1_post),
    median = apply(admin1_post, 2, stats::median),
    sd     = apply(admin1_post, 2, stats::sd),
    lower  = apply(admin1_post, 2, stats::quantile, probs = p_lo),
    upper  = apply(admin1_post, 2, stats::quantile, probs = p_hi),
    row.names = NULL
  )
  
  # Return everything you might want to reuse
  list(
    admin2_post_aligned = admin2_post,
    admin2_pop_aligned  = admin2_pop_ord,
    admin2_weights      = admin2_w,
    admin1_names        = admin1_names,
    admin1_post         = admin1_post,
    admin1_summary      = admin1_summary,
    weight_check        = check_w,
    target_order        = target_order
  )
}

## aggregate FH from admin2 to admin1
smth_agg_res1 <- aggregate_admin2_to_admin1(smth_res_ad2, admin2_pop)
smth_agg_res1_fix <- aggregate_admin2_to_admin1(smth_res_ad2_fix, admin2_pop)


## 8.3 Prepare plotting dataframe ----------
std_admin1 <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[-_/]", " ") %>%  # treat - _ / as spaces
    str_replace_all("\\s+", " ")       # squish spaces
}


df_direct <- res_ad1[["res.admin1"]] %>%
  mutate(
    admin1_std = std_admin1(admin1.name),
    region = admin1.name,
    direct = direct.est
  ) %>%
  select(admin1_std, region, direct)


df_nonfix <- smth_agg_res1$admin1_summary %>%
  mutate(
    admin1_std = std_admin1(admin1.name),
    nonfixed = mean
  ) %>%
  select(admin1_std, nonfixed)

df_fix <- smth_agg_res1_fix$admin1_summary %>%
  mutate(
    admin1_std = std_admin1(admin1.name),
    fixed = mean
  ) %>%
  select(admin1_std, fixed)

df_plot <- df_direct %>%
  left_join(df_nonfix, by = "admin1_std") %>%
  left_join(df_fix,    by = "admin1_std")


df_plot %>% filter(is.na(nonfixed) | is.na(fixed))
df_long <- df_plot %>%
  pivot_longer(cols = c("nonfixed", "fixed"),
               names_to = "type",
               values_to = "estimate")


## 8.4 Plot figure9 --------
all_vals <- c(df_plot$direct, df_plot$nonfixed, df_plot$fixed)
axis_min <- min(all_vals, na.rm = TRUE)
axis_max <- max(all_vals, na.rm = TRUE)


nested_adm1_compare <- ggplot(df_long, aes(x = direct, y = estimate, label = region)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggrepel::geom_text_repel(size = 4) +
  coord_fixed() +
  scale_x_continuous(limits = c(axis_min, axis_max)) +
  scale_y_continuous(limits = c(axis_min, axis_max)) +
  facet_wrap(~ type, nrow = 1,
             labeller = as_labeller(c(nonfixed = "Non-fixed", fixed = "Fixed"))) +
  theme_bw(base_size = 14) +
  labs(x = "Direct estimate",
       y = "Model estimate")


nested_adm1_compare

ggsave(
  filename = here::here("figures", "dhs", "nested_adm1_compare.pdf"),
  plot = nested_adm1_compare,
  width = 10,
  height = 8
)


# **********************
# 9. Table2: Admin-1 estimates table  -------
# **********************

# admin1_std: Title Case + keep hyphen (e.g., "north western" -> "North Western")
to_admin1_std <- function(x) {
  x <- gsub("[_-]+", " ", x)             # treat "_" and "-" as separators
  x <- trimws(gsub("\\s+", " ", x))      # clean extra spaces
  tools::toTitleCase(tolower(x))         # Title Case
}

ci_df <- function(post_mat, prefix) {
  apply(post_mat, 2, quantile, probs = c(0.025, 0.975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste0(c("lower_", "upper_"), prefix)) |>
    mutate(admin1_std = to_admin1_std(colnames(post_mat)))
}

df_plot_ci <- df_plot %>%
  mutate(admin1_std = to_admin1_std(admin1_std)) %>%
  left_join(ci_df(smth_agg_res1$admin1_post,     "nonfixed"), by = "admin1_std") %>%
  left_join(ci_df(smth_agg_res1_fix$admin1_post, "fixed"),    by = "admin1_std") %>%
  left_join(
    res_ad1$res.admin1 %>%
      transmute(
        admin1_std = to_admin1_std(admin1.name),
        lower_direct = direct.lower,
        upper_direct = direct.upper
      ),
    by = "admin1_std"
  ) %>%
  select(
    admin1_std,
    direct, lower_direct, upper_direct,
    nonfixed, lower_nonfixed, upper_nonfixed,
    fixed, lower_fixed, upper_fixed
  )


df_plot_ci


# **********************
# 10. National aggregate estimate -------
# **********************

std_admin1 <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[-_/]", " ") %>%
    str_squish()
}

make_qa_vec <- function(admin1_names_post, qa_table, key = "admin1_std") {
  idx <- match(std_admin1(admin1_names_post), qa_table[[key]])
  qa_vec <- qa_table$q_a[idx]
  if (anyNA(qa_vec)) {
    missing <- admin1_names_post[is.na(qa_vec)]
    stop("Missing q_a for: ", paste(missing, collapse = ", "))
  }
  qa_vec
}

nat_sum_from_posts <- function(posts, qa_vec, method_names) {
  nat_posts <- lapply(posts, \(mat) as.numeric(mat %*% qa_vec))
  
  tibble(
    method = method_names,
    mean   = sapply(nat_posts, mean),
    sd     = sapply(nat_posts, sd),
    lower  = sapply(nat_posts, \(x) as.numeric(quantile(x, 0.025))),
    upper  = sapply(nat_posts, \(x) as.numeric(quantile(x, 0.975)))
  )
}

# posterior matrices
post_direct <- res_ad1$admin1_post
post_nonfix <- smth_agg_res1$admin1_post
post_fix    <- smth_agg_res1_fix$admin1_post

admin1_names_post <- colnames(post_nonfix)  # choose one as the reference ordering

## 10.1 Survey-weighted fractions ----
qa_svy <- data0 %>%
  transmute(admin1_std = std_admin1(v024), weight) %>%
  group_by(admin1_std) %>%
  summarise(w_sum = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  mutate(q_a = w_sum / sum(w_sum)) %>%
  select(admin1_std, q_a)

qa_vec_svy <- make_qa_vec(admin1_names_post, qa_svy)

nat_sum_svy <- nat_sum_from_posts(
  posts = list(post_direct, post_nonfix, post_fix),
  qa_vec = qa_vec_svy,
  method_names = c(
    "Direct (survey-weighted)",
    "Nested unfixed (survey-weighted)",
    "Nested fixed (survey-weighted)"
  )
)

nat_sum_svy

## 10.2 WorldPop fractions ----
qa_wp <- admin1_pop %>%
  transmute(admin1_std = std_admin1(admin1.name), pop = total_pop_admin1) %>%
  mutate(q_a = pop / sum(pop)) %>%
  select(admin1_std, q_a)

qa_vec_wp <- make_qa_vec(admin1_names_post, qa_wp)

nat_sum_wp <- nat_sum_from_posts(
  posts = list(post_direct, post_nonfix, post_fix),
  qa_vec = qa_vec_wp,
  method_names = c(
    "Direct (WorldPop)",
    "Nested unfixed (WorldPop)",
    "Nested fixed (WorldPop)"
  )
)

nat_sum_wp


# **********************
# 11. Figure 10 (optional) -------
# **********************

df_diff <- df_plot %>%
  mutate(
    diff_nonfixed = nonfixed - direct,
    diff_fixed    = fixed    - direct
  )


all_vals <- c(df_diff$diff_nonfixed, df_diff$diff_fixed)
max_abs <- max(abs(all_vals))
pad <- 0.0025 
lim <- c(-max_abs - pad, max_abs + pad)

compare_adm1_diff <- ggplot(df_diff,
                            aes(x = diff_nonfixed,
                                y = diff_fixed,
                                label = region)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    size = 4,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    seed = 123
  ) +
  coord_fixed() +
  scale_x_continuous(limits = lim) +
  scale_y_continuous(limits = lim) +
  labs(
    x = "Non-fixed minus Direct",
    y = "Fixed minus Direct",
    title = "Difference from Direct Estimate: Non-fixed vs. Fixed"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10))
  )

compare_adm1_diff

ggsave(
  filename = here::here("figures", "dhs", "nested_adm1_compare_diff.pdf"),
  plot = compare_adm1_diff,
  width = 8,
  height = 8
)



