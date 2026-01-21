# *******************
# date: 2026.01.21
# task: This script is for interval plot with different models (fh is nested), 
# author: Jitong Jiang
# ********************

library(ggplot2)
library(dplyr)
library(patchwork)
library(surveyPrev)
library(SUMMER)


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

res_ad1 <- directEST_1030(data = data,
                          cluster.info = cluster.info,
                          admin = 1)

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

## nested, unfixed
smth_res_ad2 <- fhModel_1030(subset(data, !cluster %in% bad_clusters),
                        cluster.info = cluster.info,
                        admin.info = admin.info2,
                        admin = 2,
                        model = "bym2",
                        aggregation = FALSE,
                        var.fix = FALSE,
                        nested = TRUE)

## nested, fix
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

# country_shp_analysis <- readRDS("~/Desktop/vairance fix/Data/Zambia/country_shp_analysis.rds")
country_shp_analysis <- readRDS(here::here("data", "Zambia", "country_shp_analysis.rds"))

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

missing_admin2 <- unsampled_admin2$admin2.name.full


# Extract the admin2 names with low variance, good variance, and missing values.
lowVar_admin2 <- bad_admin2

goodVar_admin2 <- res_ad2$res.admin2 %>%
  filter(!(admin2.name.full %in% lowVar_admin2)) %>%
  pull(admin2.name.full)


# Using existing regional classification definitions
area_types <- data.frame(
  admin2.name.full = c(lowVar_admin2, goodVar_admin2, missing_admin2),
  area_type = c(
    rep("low variance", length(lowVar_admin2)),
    rep("usable variance", length(goodVar_admin2)),
    rep("missing", length(missing_admin2))
  )
)

# Add region type to all results
add_area_type <- function(df) {
  if(!"admin2.name.full" %in% names(df)) {
    stop("Data frame does not contain 'admin2.name.full' column")
  }
  left_join(df, area_types, by = "admin2.name.full") %>%
    mutate(area_type = factor(area_type, 
                              levels = c("legal variance", "illegal variance", "missing")))
}

# Ensure all dataframes have the 'admin2.name.full' column.
res_ad2$res.admin2 <- add_area_type(res_ad2$res.admin2)
res_ad2_fix$res.admin2 <- add_area_type(res_ad2_fix$res.admin2)
smth_res_ad2$res.admin2 <- add_area_type(smth_res_ad2$res.admin2)
smth_res_ad2_fix$res.admin2 <- add_area_type(smth_res_ad2_fix$res.admin2)


# **********************
# 5. plotting -------
# **********************


## 5.1. prepare plotting data -------

combined_data <- bind_rows(
  res_ad2$res.admin2 %>%
    select(admin1.name, admin2.name.full, area_type,
           estimate = direct.est,
           lower = direct.lower,
           upper = direct.upper) %>%
    mutate(model = "Direct"),
  
  res_ad2_fix$res.admin2 %>%
    select(admin1.name, admin2.name.full, area_type,
           estimate = direct.est,
           lower = direct.lower,
           upper = direct.upper) %>%
    mutate(model = "Fix Direct"),
  
  smth_res_ad2$res.admin2 %>%
    select(admin1.name, admin2.name.full, area_type,
           estimate = mean, lower, upper) %>%
    mutate(model = "FH"),
  
  smth_res_ad2_fix$res.admin2 %>%
    select(admin1.name, admin2.name.full, area_type,
           estimate = mean, lower, upper) %>%
    mutate(model = "Fixed FH")
) %>%
  mutate(
    area_type = case_when(
      admin2.name.full %in% lowVar_admin2 ~ "illegal variance",
      admin2.name.full %in% missing_admin2 ~ "missing",
      TRUE ~ "legal variance"
    ) %>%
      factor(levels = c("legal variance", "illegal variance", "missing")),
    
    model = factor(model, 
                   levels = c("Direct", "Fix Direct", "FH", "Fixed FH"),
                   ordered = TRUE)
  )

# get median for each Admin1
admin1_medians <- res_ad1$res.admin1 %>%
  select(admin1.name, admin1_median = direct.est)

## 5.2. plotting function -------
create_admin1_plot <- function(admin1_name) {
  plot_data <- combined_data %>%
    filter(admin1.name == admin1_name) %>%
    left_join(admin1_medians, by = "admin1.name")
  
  # Median estimate for the current Admin1
  current_median <- unique(plot_data$admin1_median)
  
  ggplot(plot_data, aes(x = reorder(admin2.name.full, as.numeric(area_type)))) +
    geom_rect(
      data = plot_data %>%
        distinct(area_type) %>%
        mutate(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = area_type),
      inherit.aes = FALSE,
      alpha = 0.15
    ) +
    
    # Reference line for the Admin1 median
    geom_hline(
      yintercept = current_median,
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    
    # Point estimates and 95% CIs for each model
    geom_pointrange(
      aes(y = estimate, ymin = lower, ymax = upper,
          color = model, shape = model),
      position = position_dodge(width = 0.8),
      size = 0.5,
      fatten = 2
    ) +
    
    # Color mapping
    scale_color_manual(
      values = c(
        "Direct"    = "#ff7f0e",
        "Fix Direct" = "#f768a1",
        "FH"        = "#9467bd",
        "Fixed FH"  = "#2ca02c"
      )
    ) +
    
    # Shape mapping
    scale_shape_manual(
      values = c(
        "Direct"    = 17,
        "Fix Direct" = 18,
        "FH"        = 3,
        "Fixed FH"  = 15
      )
    ) +
    
    # Background fill for area types (no legend)
    scale_fill_manual(
      values = c(
        "usable variance" = "white",
        "low variance"    = "lightgrey",
        "missing"         = "grey70"
      ),
      guide = "none"
    ) +
    
    # Labels and theme settings
    labs(
      title    = paste("Admin1:", admin1_name),
      subtitle = paste("Admin1 Median:", round(current_median, 4)),
      x        = NULL,
      y        = "Estimate with 95% CI",
      color    = "Model",
      shape    = "Model"
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position    = "bottom",
      plot.title         = element_text(size = 12, face = "bold"),
      plot.subtitle      = element_text(color = "red", size = 10),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold")
    ) +
    
    # Facet by area_type
    facet_grid(. ~ area_type, scales = "free_x", space = "free_x")
}


## 5.3. plots for each admin1 -------

admin_list <- unique(combined_data$admin1.name)

for(admin in admin_list) {
  p <- create_admin1_plot(admin)
  print(p)
  ggsave(
    filename = here::here("figures", "dhs",
      paste0("Nested_Admin1_", admin, ".png")
    ),
    plot = p,
    width = 8,     
    height = 6,    
    dpi = 300,
    bg = "white"
  )
}


## 5.4. combine all the plots -------

plot_list <- lapply(admin_list, create_admin1_plot)

combined_plot <- wrap_plots(plot_list, ncol = 2, nrow = 5) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

combined_plot

ggsave(
  filename = here::here(
    "figures", "dhs", "Nested_Admin1_all.png"
  ),
  plot = combined_plot,
  width = 14,   
  height = 18,  
  dpi = 300,
  bg = "white"
)

