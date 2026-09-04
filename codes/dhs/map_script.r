# *******************
# date: 2026.01.21
# task: This script is for map plots for direct model and FH model
# author: Jitong Jiang
# ********************

library(surveyPrev)
library(SUMMER)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(grid)
library(here)


source(here::here("codes", "directEST_1030_national.R"))
source(here::here("codes", "fhModel_1030.R"))


# **********************
# 1. Basic setting ----
# **********************

indicator <- "CN_NUTS_C_WH2"

year <- 2018
frame_year <- 2010
country <- "Zambia"
country.abbrev <- "ZMB"


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
  filter(!is.na(value))


geo <- getDHSgeo(
  country = country,
  year = year
)

poly.adm1 <- geodata::gadm(
  country = country.abbrev,
  level = 1,
  path = tempdir()
)

poly.adm1 <- sf::st_as_sf(poly.adm1)


poly.adm2 <- geodata::gadm(
  country = country.abbrev,
  level = 2,
  path = tempdir()
)

poly.adm2 <- sf::st_as_sf(poly.adm2) %>%
  mutate(
    admin2.name.full = paste0(NAME_1, "_", NAME_2)
  )


cluster.info <- clusterInfo(
  geo = geo,
  poly.adm1 = poly.adm1,
  poly.adm2 = poly.adm2,
  by.adm1 = "NAME_1",
  by.adm2 = "NAME_2"
)

admin.info1 <- adminInfo(
  poly.adm = poly.adm1,
  admin = 1,
  by.adm = "NAME_1"
)

admin.info2 <- adminInfo(
  poly.adm = poly.adm2,
  admin = 2,
  by.adm = "NAME_2",
  by.adm.upper = "NAME_1"
)


# **********************
# 2. Direct estimate results ----
# **********************

res_ad2 <- directEST_1030(
  data = data,
  cluster.info = cluster.info,
  admin = 2,
  aggregation = FALSE,
  var.fix = FALSE
)


res_ad2_fix <- directEST_1030(
  data = data,
  cluster.info = cluster.info,
  admin = 2,
  aggregation = FALSE,
  var.fix = TRUE
)


# **********************
# 3. FH estimate results ----
# **********************

bad_admin2 <- res_ad2_fix$fixed_areas

bad_clusters <- subset(
  cluster.info$data,
  admin2.name.full %in% bad_admin2
)$cluster


smth_res_ad2 <- fhModel_1030(
  subset(data, !cluster %in% bad_clusters),
  cluster.info = cluster.info,
  admin.info = admin.info2,
  admin = 2,
  model = "bym2",
  aggregation = FALSE,
  var.fix = FALSE,
  nested = TRUE
)


smth_res_ad2_fix <- fhModel_1030(
  data,
  cluster.info = cluster.info,
  admin.info = admin.info2,
  admin = 2,
  model = "bym2",
  aggregation = FALSE,
  var.fix = TRUE,
  nested = TRUE
)


# **********************
# 4. Classification by type ----
# **********************

country_shp_analysis <- readRDS(
  "~/Desktop/vairance fix/Data/Zambia/country_shp_analysis.rds"
)


data0 <- data %>%
  filter(!is.na(value))


myData_tmp <- data0 %>%
  group_by(cluster) %>%
  dplyr::summarise(
    strata = unique(strata),
    Ntrials = n(),
    value = sum(value),
    households_number = n_distinct(householdID),
    .groups = "drop"
  )


myData <- merge(
  myData_tmp,
  cluster.info[["data"]],
  by = "cluster"
)


sampled_admin2_cluster <- myData %>%
  group_by(
    admin1.name,
    admin2.name,
    admin2.name.full
  ) %>%
  summarise(
    sampled_admin2_cluster = n(),
    .groups = "drop"
  )


full_admin2 <- country_shp_analysis[["Admin-2"]] |>
  dplyr::mutate(
    admin2.name.full = paste(
      NAME_1,
      NAME_2,
      sep = "_"
    )
  )


unsampled_admin2 <- dplyr::filter(
  full_admin2,
  !(admin2.name.full %in%
      sampled_admin2_cluster$admin2.name.full)
)


sparse_admin2 <- bad_admin2

missing_admin2 <- unsampled_admin2$admin2.name.full


# **********************
# 5. Prepare mapping datasets ----
# **********************

res_map_df <- poly.adm2 %>%
  left_join(
    res_ad2[["res.admin2"]],
    by = "admin2.name.full"
  ) %>%
  mutate(
    CI_width = direct.upper - direct.lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )


res_map_df_fix <- poly.adm2 %>%
  left_join(
    res_ad2_fix[["res.admin2"]],
    by = "admin2.name.full"
  ) %>%
  mutate(
    CI_width = direct.upper - direct.lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )


smth_map_df <- poly.adm2 %>%
  left_join(
    smth_res_ad2[["res.admin2"]],
    by = "admin2.name.full"
  ) %>%
  mutate(
    CI_width = upper - lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )


smth_fix_map_df <- poly.adm2 %>%
  left_join(
    smth_res_ad2_fix[["res.admin2"]],
    by = "admin2.name.full"
  ) %>%
  mutate(
    CI_width = upper - lower,
    hatch_type = case_when(
      admin2.name.full %in% missing_admin2 ~ "missing",
      admin2.name.full %in% sparse_admin2 ~ "sparse",
      TRUE ~ "none"
    )
  )


# **********************
# 6. Common color scales ----
# **********************

## 6.1 Shared range for point estimates

all_point_estimates <- c(
  res_map_df$direct.est,
  res_map_df_fix$direct.est,
  smth_map_df$mean,
  smth_fix_map_df$mean
)

point_limits <- range(
  all_point_estimates,
  na.rm = TRUE
)


## 6.2 Shared range for CI widths

all_ci_widths <- c(
  res_map_df$CI_width,
  res_map_df_fix$CI_width,
  smth_map_df$CI_width,
  smth_fix_map_df$CI_width
)

ci_limits <- range(
  all_ci_widths,
  na.rm = TRUE
)


# Illegal areas

illegal_areas <- res_map_df %>%
  filter(hatch_type == "sparse") %>%
  mutate(
    area_number = seq_len(n())
  )


# **********************
# 7. Mapping function ----
# **********************

create_map_with_legend <- function(
    data,
    fill_var,
    fill_name,
    fill_limits,
    show_legend = TRUE,
    legend_position = "bottom",
    illegal_areas = NULL) {
  
  if (is.null(illegal_areas)) {
    
    illegal_areas <- data %>%
      filter(hatch_type == "sparse") %>%
      mutate(
        area_number = seq_len(n())
      )
  }
  
  
  p <- ggplot(data) +
    
    geom_sf(
      aes(fill = {{ fill_var }}),
      color = "black",
      size = 0.1,
      na.rm = TRUE
    ) +
    
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
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      )
    )
  
  
  if (show_legend) {
    
    p <- p +
      theme(
        legend.position = legend_position,
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0)
      )
    
  } else {
    
    p <- p +
      theme(
        legend.position = "none"
      )
  }
  
  
  return(p)
}


# **********************
# 8. Create 8 map panels ----
#    All legends removed
# **********************

## 8.1 Direct ----

p_res_mean_panel <- create_map_with_legend(
  data = res_map_df,
  fill_var = direct.est,
  fill_name = "Estimate",
  fill_limits = point_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


p_res_ci_panel <- create_map_with_legend(
  data = res_map_df,
  fill_var = CI_width,
  fill_name = "CI Width",
  fill_limits = ci_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


## 8.2 Fixed direct ----

p_res_fix_mean_panel <- create_map_with_legend(
  data = res_map_df_fix,
  fill_var = direct.est,
  fill_name = "Estimate",
  fill_limits = point_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


p_res_fix_ci_panel <- create_map_with_legend(
  data = res_map_df_fix,
  fill_var = CI_width,
  fill_name = "CI Width",
  fill_limits = ci_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


## 8.3 Fay-Herriot ----

p_fh_mean_panel <- create_map_with_legend(
  data = smth_map_df,
  fill_var = mean,
  fill_name = "Estimate",
  fill_limits = point_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


p_fh_ci_panel <- create_map_with_legend(
  data = smth_map_df,
  fill_var = CI_width,
  fill_name = "CI Width",
  fill_limits = ci_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


## 8.4 Fixed Fay-Herriot ----

p_fh_fix_mean_panel <- create_map_with_legend(
  data = smth_fix_map_df,
  fill_var = mean,
  fill_name = "Estimate",
  fill_limits = point_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


p_fh_fix_ci_panel <- create_map_with_legend(
  data = smth_fix_map_df,
  fill_var = CI_width,
  fill_name = "CI Width",
  fill_limits = ci_limits,
  show_legend = FALSE,
  illegal_areas = illegal_areas
) +
  labs(title = NULL)


# **********************
# 9. Separate color bars ----
# **********************

## 9.1 Estimate horizontal color bar ----

legend_mean_horizontal <- get_legend(
  
  create_map_with_legend(
    data = smth_fix_map_df,
    fill_var = mean,
    fill_name = "Estimate",
    fill_limits = point_limits,
    show_legend = TRUE,
    legend_position = "bottom",
    illegal_areas = illegal_areas
  ) +
    
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(4, "cm"),
        barheight = unit(0.35, "cm")
      )
    ) +
    
    theme(
      legend.position = "bottom"
    )
)


## 9.2 Estimate vertical color bar ----

legend_mean_vertical <- get_legend(
  
  create_map_with_legend(
    data = smth_fix_map_df,
    fill_var = mean,
    fill_name = "Estimate",
    fill_limits = point_limits,
    show_legend = TRUE,
    legend_position = "right",
    illegal_areas = illegal_areas
  ) +
    
    guides(
      fill = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(0.35, "cm"),
        barheight = unit(4, "cm")
      )
    ) +
    
    theme(
      legend.position = "right"
    )
)


## 9.3 CI Width horizontal color bar -----

legend_ci_horizontal <- get_legend(
  
  create_map_with_legend(
    data = smth_fix_map_df,
    fill_var = CI_width,
    fill_name = "CI Width",
    fill_limits = ci_limits,
    show_legend = TRUE,
    legend_position = "bottom",
    illegal_areas = illegal_areas
  ) +
    
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(4, "cm"),
        barheight = unit(0.35, "cm")
      )
    ) +
    
    theme(
      legend.position = "bottom"
    )
)


## 9.4 CI Width vertical color bar -----

legend_ci_vertical <- get_legend(
  
  create_map_with_legend(
    data = smth_fix_map_df,
    fill_var = CI_width,
    fill_name = "CI Width",
    fill_limits = ci_limits,
    show_legend = TRUE,
    legend_position = "right",
    illegal_areas = illegal_areas
  ) +
    
    guides(
      fill = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(0.35, "cm"),
        barheight = unit(4, "cm")
      )
    ) +
    
    theme(
      legend.position = "right"
    )
)


# **********************
# 10. Convert legends to plots -------
# **********************

legend_mean_horizontal_plot <- ggdraw() +
  draw_grob(legend_mean_horizontal)


legend_mean_vertical_plot <- ggdraw() +
  draw_grob(legend_mean_vertical)


legend_ci_horizontal_plot <- ggdraw() +
  draw_grob(legend_ci_horizontal)


legend_ci_vertical_plot <- ggdraw() +
  draw_grob(legend_ci_vertical)


# **********************
# 11. Save color bars -----
# **********************

ggsave(
  filename = here::here(
    "figures", "dhs", "new",
    "colorbar_estimate_horizontal.png"
  ),
  plot = legend_mean_horizontal_plot,
  width = 5,
  height = 1,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = here::here(
    "figures", "dhs", "new",
    "colorbar_estimate_vertical.png"
  ),
  plot = legend_mean_vertical_plot,
  width = 1.5,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = here::here(
    "figures", "dhs", "new",
    "colorbar_CI_width_horizontal.png"
  ),
  plot = legend_ci_horizontal_plot,
  width = 5,
  height = 1,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = here::here(
    "figures", "dhs", "new",
    "colorbar_CI_width_vertical.png"
  ),
  plot = legend_ci_vertical_plot,
  width = 1.5,
  height = 5,
  dpi = 300,
  bg = "white"
)

# **********************
# 12. Optional: combined 8-panel plot without legends ------
# **********************

row_lab <- function(label) {
  
  ggdraw() +
    
    draw_label(
      label,
      fontface = "bold",
      x = 0.5,
      y = 0.5,
      angle = 90,
      size = 10
    ) +
    
    theme(
      plot.margin = margin(
        0,
        0.1,
        0,
        0.1
      )
    )
}


row_direct <- plot_grid(
  row_lab("Direct"),
  p_res_mean_panel,
  p_res_ci_panel,
  nrow = 1,
  rel_widths = c(
    0.05,
    1,
    1
  )
)


row_direct_fix <- plot_grid(
  row_lab("Fixed direct"),
  p_res_fix_mean_panel,
  p_res_fix_ci_panel,
  nrow = 1,
  rel_widths = c(
    0.05,
    1,
    1
  )
)


row_fh <- plot_grid(
  row_lab("Fay-Herriot"),
  p_fh_mean_panel,
  p_fh_ci_panel,
  nrow = 1,
  rel_widths = c(
    0.05,
    1,
    1
  )
)


row_fh_fix <- plot_grid(
  row_lab("Fixed Fay-Herriot"),
  p_fh_fix_mean_panel,
  p_fh_fix_ci_panel,
  nrow = 1,
  rel_widths = c(
    0.05,
    1,
    1
  )
)


title_row <- plot_grid(
  ggdraw(),
  ggdraw() +
    draw_label(
      "Mean",
      fontface = "bold",
      size = 11
    ),
  ggdraw() +
    draw_label(
      "CI Width",
      fontface = "bold",
      size = 11
    ),
  nrow = 1,
  rel_widths = c(
    0.05,
    1,
    1
  )
)


combined_plot_4x2 <- plot_grid(
  title_row,
  row_direct,
  row_direct_fix,
  row_fh,
  row_fh_fix,
  ncol = 1,
  rel_heights = c(
    0.08,
    1,
    1,
    1,
    1
  )
)


print(combined_plot_4x2)


# **********************
# Save 8 individual maps -----
# **********************

ggsave(
  here::here("figures", "dhs", "new", "Direct_Estimate.png"),
  p_res_mean_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

ggsave(
  here::here("figures", "dhs", "new", "Direct_CI_Width.png"),
  p_res_ci_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

ggsave(
  here::here("figures", "dhs", "new", "Fixed_Direct_Estimate.png"),
  p_res_fix_mean_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

ggsave(
  here::here("figures", "dhs", "new", "Fixed_Direct_CI_Width.png"),
  p_res_fix_ci_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

ggsave(
  here::here("figures", "dhs", "new", "FH_Estimate.png"),
  p_fh_mean_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

ggsave(
  here::here("figures", "dhs", "new", "FH_CI_Width.png"),
  p_fh_ci_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

ggsave(
  here::here("figures", "dhs", "new", "Fixed_FH_Estimate.png"),
  p_fh_fix_mean_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

ggsave(
  here::here("figures", "dhs", "new", "Fixed_FH_CI_Width.png"),
  p_fh_fix_ci_panel,
  width = 5, height = 5, dpi = 300, bg = "white"
)

# **********************
# 13. Optional checks ------
# **********************

print(legend_mean_horizontal_plot)
print(legend_mean_vertical_plot)
print(legend_ci_horizontal_plot)
print(legend_ci_vertical_plot)


