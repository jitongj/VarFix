# *******************
# date: 2026.01.21
# task: This script is for map plots for direct model and FH model, 
# author: Jitong Jiang
# ********************

library(surveyPrev)
library(SUMMER)
library(dplyr)
library(ggplot2)
library(cowplot)  
library(viridis)
library(grid) 


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
# poly.adm2=poly.adm2[poly.adm2$ENGTYPE_2=="Local Authority",]


cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, 
                            by.adm1 = "NAME_1",by.adm2 = "NAME_2")
admin.info1 <- adminInfo(poly.adm = poly.adm1, admin = 1, by.adm = "NAME_1")
admin.info2 <- adminInfo(poly.adm = poly.adm2, admin = 2,
                         by.adm = "NAME_2", by.adm.upper = "NAME_1")


# **********************
# 2. Direct estimate results-------
# **********************
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


smth_res_ad2_fix <- fhModel_1030(data,
                            cluster.info = cluster.info,
                            admin.info = admin.info2,
                            admin = 2,
                            model = "bym2",
                            aggregation = FALSE,
                            var.fix = TRUE,
                            nested = TRUE)

# **********************
# 4. Classification by type -------
# **********************

country_shp_analysis <- readRDS("~/Desktop/vairance fix/Data/Zambia/country_shp_analysis.rds")

data0 <- data %>% filter(!is.na(value))

myData_tmp <- data0 %>%
  group_by(cluster) %>%   # Group by cluster
  dplyr::summarise(       # Compute within-cluster summary statistics
    strata = unique(strata),  # Unique strata value for the cluster (assuming one strata per cluster)
    Ntrials = n(),            # Number of individuals in the cluster
    value = sum(value),       # Total number of "1"s (sum of individual values) in the cluster
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
# 5. prepare mapping dataset -------
# **********************

res_map_df <- poly.adm2 %>%
  left_join(res_ad2[["res.admin2"]], by = "admin2.name.full") %>%
  mutate(
    CI_width = direct.upper - direct.lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )

res_map_df_fix <- poly.adm2 %>%
  left_join(res_ad2_fix[["res.admin2"]], by = "admin2.name.full") %>%
  mutate(
    CI_width = direct.upper - direct.lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )

smth_map_df <- poly.adm2 %>%
  left_join(smth_res_ad2[["res.admin2"]], by = "admin2.name.full") %>%
  mutate(
    CI_width = upper - lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )

smth_fix_map_df <- poly.adm2 %>%
  left_join(smth_res_ad2_fix[["res.admin2"]], by = "admin2.name.full") %>%
  mutate(
    CI_width = upper - lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )


# Compute common color scale limits across all models
# 1. Shared range for all point estimates
all_point_estimates <- c(
  res_map_df$direct.est, 
  res_map_df_fix$direct.est,
  smth_map_df$mean, 
  smth_fix_map_df$mean
)
point_limits <- range(all_point_estimates, na.rm = TRUE)

# 2. Shared range for all CI widths
all_ci_widths <- c(
  res_map_df$CI_width, 
  res_map_df_fix$CI_width,
  smth_map_df$CI_width, 
  smth_fix_map_df$CI_width
)
ci_limits <- range(all_ci_widths, na.rm = TRUE)


# Extract illegal_areas information from the original results
illegal_areas <- res_map_df %>%
  filter(hatch_type == "sparse") %>%
  mutate(area_number = seq_len(n()))

# **********************
# 6. mapping function -------
# **********************
create_map_with_legend <- function(data, fill_var, fill_name, fill_limits, 
                                   show_legend = TRUE, legend_position = "bottom",
                                   illegal_areas = NULL) {
  
  if (is.null(illegal_areas)) {
    illegal_areas <- data %>%
      filter(hatch_type == "sparse") %>%
      mutate(area_number = seq_len(n()))
  }
  
  p <- ggplot(data) +
    geom_sf(aes(fill = {{fill_var}}), color = "black", size = 0.1, na.rm = TRUE) +
    geom_sf(
      data = illegal_areas,
      fill = NA, 
      color = "red", 
      size = 0.8,
      linetype = "solid"
    ) +
    geom_sf_text(
      data = illegal_areas,
      aes(label = area_number),
      color = "red",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      direction = -1,
      name = fill_name,
      limits = fill_limits,
      na.value = "grey90"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7),
      plot.margin = margin(2, 2, 2, 2),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  if (show_legend) {
    p <- p +
      theme(
        legend.position = legend_position,
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0)
      )
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}


# **********************
# 7. plot maps -------
# **********************

## 7.1 Create all map panels ----

# Mean map (no legend)
p_res_mean_panel <- create_map_with_legend(
  data = res_map_df, fill_var = direct.est, fill_name = "Estimate",
  fill_limits = point_limits, show_legend = FALSE,
  illegal_areas = illegal_areas
) + labs(title = NULL)

p_res_fix_mean_panel <- create_map_with_legend(
  data = res_map_df_fix, fill_var = direct.est, fill_name = "Estimate", 
  fill_limits = point_limits, show_legend = FALSE,
  illegal_areas = illegal_areas
) + labs(title = NULL)

p_fh_mean_panel <- create_map_with_legend(
  data = smth_map_df, fill_var = mean, fill_name = "Estimate",
  fill_limits = point_limits, show_legend = FALSE,
  illegal_areas = illegal_areas
) + labs(title = NULL)

# CI width map (no legend)
p_res_ci_panel <- create_map_with_legend(
  data = res_map_df, fill_var = CI_width, fill_name = "CI Width",
  fill_limits = ci_limits, show_legend = FALSE,
  illegal_areas = illegal_areas
) + labs(title = NULL)

p_res_fix_ci_panel <- create_map_with_legend(
  data = res_map_df_fix, fill_var = CI_width, fill_name = "CI Width",
  fill_limits = ci_limits, show_legend = FALSE,
  illegal_areas = illegal_areas
) + labs(title = NULL)

p_fh_ci_panel <- create_map_with_legend(
  data = smth_map_df, fill_var = CI_width, fill_name = "CI Width",
  fill_limits = ci_limits, show_legend = FALSE,
  illegal_areas = illegal_areas
) + labs(title = NULL)


# p_fh_fix_mean_with_legend <- create_map_with_legend(
#   data = smth_fix_map_df, fill_var = mean, fill_name = "Estimate",
#   fill_limits = point_limits, show_legend = TRUE,
#   illegal_areas = illegal_areas
# ) + labs(title = NULL)
# 
# p_fh_fix_ci_with_legend <- create_map_with_legend(
#   data = smth_fix_map_df, fill_var = CI_width, fill_name = "CI Width", 
#   fill_limits = ci_limits, show_legend = TRUE,
#   illegal_areas = illegal_areas
# ) + labs(title = NULL)

p_fh_fix_mean_panel <- create_map_with_legend(
  data = smth_fix_map_df, 
  fill_var = mean, 
  fill_name = "Estimate",
  fill_limits = point_limits, 
  show_legend = FALSE,
  illegal_areas = illegal_areas
)

p_fh_fix_ci_panel <- create_map_with_legend(
  data = smth_fix_map_df, 
  fill_var = CI_width, 
  fill_name = "CI Width", 
  fill_limits = ci_limits, 
  show_legend = FALSE,
  illegal_areas = illegal_areas
)


legend_mean <- get_legend(
  create_map_with_legend(
    data = smth_fix_map_df, 
    fill_var = mean, 
    fill_name = "Estimate",
    fill_limits = point_limits, 
    show_legend = TRUE
  ) + theme(legend.position = "bottom")
)

legend_ci <- get_legend(
  create_map_with_legend(
    data = smth_fix_map_df, 
    fill_var = CI_width, 
    fill_name = "CI Width",
    fill_limits = ci_limits, 
    show_legend = TRUE
  ) + theme(legend.position = "bottom")
)


## 7.2 Combined all maps ----

# Row label helper function
row_lab <- function(label) {
  ggdraw() +
    draw_label(
      label,
      fontface = "bold",
      x = 0.5, y = 0.5,
      angle = 90,
      size = 10
    ) +
    theme(
      plot.margin = margin(0, 0.1, 0, 0.1)
    )
}

# Rows 1-3: no colorbars
row_direct <- plot_grid(
  row_lab("Direct"),
  p_res_mean_panel,
  p_res_ci_panel,
  nrow = 1,
  rel_widths = c(0.05, 1, 1)
)

row_direct_fix <- plot_grid(
  row_lab("Fixed direct"),
  p_res_fix_mean_panel,
  p_res_fix_ci_panel, 
  nrow = 1,
  rel_widths = c(0.05, 1, 1)
)

row_fh <- plot_grid(
  row_lab("Fay-Herriot"),
  p_fh_mean_panel,
  p_fh_ci_panel,
  nrow = 1, 
  rel_widths = c(0.05, 1, 1)
)

# 4th row: with colorbar
row_fh_fix <- plot_grid(
  row_lab("Fixed Fay-Herriot"),
  p_fh_fix_mean_panel,
  p_fh_fix_ci_panel,
  nrow = 1,
  rel_widths = c(0.05, 1, 1)
)

legend_row <- plot_grid(
  legend_mean, 
  legend_ci,
  nrow = 1,
  rel_widths = c(1, 1)
)


title_row <- plot_grid(
  ggdraw(),
  ggdraw() + draw_label("Mean", fontface = "bold", size = 11),
  ggdraw() + draw_label("CI Width", fontface = "bold", size = 11),
  nrow = 1,
  rel_widths = c(0.05, 1, 1)
)


combined_plot_4x2 <- plot_grid(
  title_row,
  row_direct,
  row_direct_fix,
  row_fh,
  row_fh_fix,      
  legend_row,     
  ncol = 1,
  rel_heights = c(0.08, 1, 1, 1, 1, 0.15)
)

print(combined_plot_4x2)
