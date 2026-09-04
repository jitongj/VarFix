# *******************
# date: 2026.01.21
# task: This script is for ranking_map plots for dhs data, 
# author: Jitong Jiang
# ********************


library(ggplot2)
library(dplyr)
library(patchwork)
library(surveyPrev)
library(SUMMER)
library(tidyr)

source(here::here("codes", "directEST_1030_national.R"))
source(here::here("codes", "fhModel_1030.R"))

## This is the dataset for plotting
load(here::here("data/ranking_map.RData"))


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
set.seed(2024)
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
# 4. ranking function-------
# **********************

cateRank <- function(model, 
                     admin = 2,
                     method = c("direct", "mbe"),
                     high_color   = "darkred",
                     middle_color = "darkorange",
                     low_color    = "darkgreen",
                     title = NULL,
                     lower = 0.20,        # bottom 20%
                     upper = 0.80,        # top 20% (1 - upper)
                     mid_keep_top = 5,
                     mid_keep_bottom = 5,
                     bad_admin2 = NULL,
                     prob_threshold = 0.5
) {
  
  method <- match.arg(method)
  
  admin_name <- paste0("admin", admin)          # "admin2"
  admin_res  <- paste0("res.", admin_name)      # "res.admin2"
  if (admin == 1) {
    admin_region_name <- paste0(admin_name, ".name")       # "admin1.name"
  } else {
    admin_region_name <- paste0(admin_name, ".name.full")  # "admin2.name.full"
  }
  
  # posterior
  posterior_name <- paste0(admin_name, "_post") # "admin2_post"
  posterior_data <- model[[posterior_name]]
  posterior_data <- as.matrix(posterior_data)
  n_draws <- nrow(posterior_data)
  n_areas <- ncol(posterior_data)
  
  area_names <- NULL
  

  if (!is.null(model[[admin_res]]) &&
      !is.null(model[[admin_res]][[admin_region_name]]) &&
      length(model[[admin_res]][[admin_region_name]]) == n_areas) {
    area_names <- as.character(model[[admin_res]][[admin_region_name]])
  }

  if (is.null(area_names) || length(area_names) != n_areas) {
    attr_name <- paste0("admin", admin, ".name.full")
    att <- attr(model, attr_name, exact = TRUE)
    if (!is.null(att) && length(att) == n_areas) {
      area_names <- as.character(att)
    }
  }
  
  if (is.null(area_names) || length(area_names) != n_areas) {
    message("Warning: admin names missing or mismatched, using generic Area_1, ..., Area_n.")
    area_names <- paste0("Area_", seq_len(n_areas))
  }
  
  area_names <- gsub("_", " - ", area_names)
  stopifnot(length(area_names) == n_areas)
  

  if (!is.null(bad_admin2)) {
    bad_admin2_clean <- gsub("_", " - ", bad_admin2)
  } else {
    bad_admin2_clean <- character(0)
  }
  

  ranks <- t(apply(posterior_data, 1, function(row)
    rank(-row, ties.method = "first")))  # 1 = max
  
  # top / bottom  proportion(default 20% / 20%)
  frac_top    <- (1 - upper)
  frac_bottom <- lower
  base_frac   <- min(frac_top, frac_bottom)
  group_k     <- max(1L, round(base_frac * n_areas))
  if (2 * group_k > n_areas) {
    group_k <- floor(n_areas / 2)
  }
  top_k    <- group_k
  bottom_k <- group_k
  
  mid_lo <- if (top_k > 0) top_k + 1L else 1L
  mid_hi <- if (bottom_k > 0) n_areas - bottom_k else n_areas
  bot_lo <- if (bottom_k > 0) mid_hi + 1L else n_areas + 1L
  
  # bin_idx: 1=Top, 2=Middle, 3=Bottom
  bin_idx <- matrix(2L, nrow = n_draws, ncol = n_areas)   # default is middle
  if (top_k > 0)    bin_idx[ranks <= top_k]   <- 1L
  if (bottom_k > 0) bin_idx[ranks >= bot_lo]  <- 3L
  
  
  Highest_25 <- colMeans(bin_idx == 1L) * 100  # P(Top frac)
  Middle_50  <- colMeans(bin_idx == 2L) * 100  # P(Mid frac)
  Lowest_25  <- colMeans(bin_idx == 3L) * 100  # P(Bottom frac)
  
  # point estimate
  if (method == "direct") {
    est <- "direct.est"
    var <- "direct.var"
  } else {
    est <- "mean"
    var <- "var"
  }
  
  all_values <- model[[admin_res]][[est]]
  prevalence <- as.numeric(all_values) * 100   # percentage
  
  summary_df_full <- data.frame(
    Area        = area_names,
    Lowest_25   = Lowest_25,
    Middle_50   = Middle_50,
    Highest_25  = Highest_25,
    prevalence  = prevalence,
    stringsAsFactors = FALSE
  )
  
  
  summary_df_full$color_group <- middle_color
  thr_pct <- prob_threshold * 100  
  

  summary_df_full$color_group[summary_df_full$Highest_25 >= thr_pct] <- high_color
  

  summary_df_full$color_group[
    summary_df_full$Lowest_25 >= thr_pct &
      summary_df_full$color_group != high_color
  ] <- low_color
  
  # ranking (high -> middel -> low)
  df_high <- dplyr::filter(summary_df_full, color_group == high_color)   |> 
    dplyr::arrange(dplyr::desc(Highest_25))
  df_mid  <- dplyr::filter(summary_df_full, color_group == middle_color) |> 
    dplyr::arrange(dplyr::desc(Highest_25), Lowest_25)
  df_low  <- dplyr::filter(summary_df_full, color_group == low_color)    |> 
    dplyr::arrange(Lowest_25)
  
  summary_df_full <- rbind(df_high, df_mid, df_low)
  
  
  n_mid <- nrow(df_mid)
  
  if (n_mid > (mid_keep_top + mid_keep_bottom)) {
    keep_top <- min(mid_keep_top, n_mid)
    keep_bot <- min(mid_keep_bottom, n_mid - keep_top)
    
    mid_top <- df_mid[1:keep_top, , drop = FALSE]
    mid_bot <- df_mid[(n_mid - keep_bot + 1):n_mid, , drop = FALSE]
    
    placeholder_row <- data.frame(
      Area        = "...",
      Lowest_25   = NA_real_,
      Middle_50   = NA_real_,
      Highest_25  = NA_real_,
      prevalence  = NA_real_,
      color_group = middle_color,
      stringsAsFactors = FALSE
    )
    
    display_df <- rbind(
      df_high,
      mid_top,
      placeholder_row,
      mid_bot,
      df_low
    )
  } else {
    display_df <- rbind(df_high, df_mid, df_low)
  }
  
  
  plot_data <- tidyr::pivot_longer(
    display_df,
    cols = c("Highest_25", "Middle_50", "Lowest_25"),
    names_to = "Category",
    values_to = "Percentage"
  )
  
  
  y_levels <- rev(display_df$Area)
  plot_data$Area     <- factor(plot_data$Area, levels = y_levels)
  plot_data$Category <- factor(
    plot_data$Category,
    levels = c("Highest_25", "Middle_50", "Lowest_25")
  )
  
  # label bad_admin2 with * in y axis
  y_labels <- ifelse(
    y_levels %in% bad_admin2_clean & y_levels != "...",
    paste0("* ", y_levels),
    y_levels
  )
  

  color_vec <- sapply(y_levels, function(a) {
    display_df$color_group[match(a, display_df$Area)]
  })
  

  n_high_theoretical <- top_k     
  n_low_theoretical <- bottom_k  
  n_mid_theoretical <- n_areas - n_high_theoretical - n_low_theoretical  
  

  n_high_actual <- sum(summary_df_full$color_group == high_color)
  n_low_actual <- sum(summary_df_full$color_group == low_color)
  n_mid_actual <- sum(summary_df_full$color_group == middle_color)
  
  x_labs <- c(
    paste0("Highest\n", round(frac_top * 100), "% (", n_high_theoretical, ")"),
    paste0("Middle\n",  round((upper - lower) * 100), "% (", n_mid_theoretical,  ")"),
    paste0("Lowest\n",  round(frac_bottom * 100), "% (", n_low_theoretical, ")")
  )
  

  p <- ggplot2::ggplot(
    data = subset(plot_data, !is.na(Percentage)),
    ggplot2::aes(x = Category, y = Area)
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill = Percentage), color = "grey", lwd = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", Percentage)),
      color = "black", size = 4
    ) +
    ggplot2::scale_fill_gradient(low = "white", high = "royalblue",
                                 name = "Percentage(%)", limits = c(0, 100)) +
    ggplot2::scale_x_discrete(labels = x_labs) +
    ggplot2::scale_y_discrete(drop = FALSE, labels = y_labels) +
    ggplot2::labs(
      title = title,
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 20
      ),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = color_vec, size = 12),
      axis.text.x = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 16),
      legend.text  = ggplot2::element_text(size = 15)
    )

  
  p <- p +
    ggplot2::geom_text(
      data = subset(display_df, !is.na(prevalence)),
      ggplot2::aes(
        x = 0,
        y = factor(Area, levels = y_levels),
        label = sprintf("%.1f", prevalence)
      ),
      hjust = 0, size = 4, inherit.aes = FALSE
    )
  

  colnames(summary_df_full) <- c(
    admin_region_name, "Lowest", "Middle", "Highest", "Prevalence", "color_group"
  )
  

  message(paste("Actual highlighted areas - High:", n_high_actual, "Low:", n_low_actual))
  
  list(
    plot  = p,
    table = summary_df_full,
    actual_counts = c(high = n_high_actual, middle = n_mid_actual, low = n_low_actual)
  )
}



# **********************
# 5. ranking plots-------
# **********************

rank_smth2 <- cateRank(
  model  = smth_res_ad2,
  admin  = 2,
  method = "mbe",
  upper  = 0.8,
  lower  = 0.2,
  title  = "Fay-Herriot Model",
  bad_admin2 = bad_admin2
)

rank_smth2$plot


rank_smth2_fix <- cateRank(
  model  = smth_res_ad2_fix,
  admin  = 2,
  method = "mbe",
  upper  = 0.8,
  lower  = 0.2,
  title  = "Fixed Fay-Herriot Model",
  bad_admin2 = bad_admin2
)

rank_smth2_fix$plot


# **********************
# 6. map function-------
# **********************

high_color   = "darkred"
middle_color = "darkorange"
low_color    = "darkgreen"

a = rank_smth2
library(dplyr)
library(ggplot2)
library(sf)


rank_map <- function(rank_obj,
                     shapefile,
                     admin_col   = "admin2.name.full",
                     high_color   = "darkred",
                     middle_color = "darkorange",
                     low_color    = "darkgreen",
                     title        = "") {
  
 
  shapefile <- shapefile %>%
    mutate(
      !!admin_col := gsub("_", " - ", .data[[admin_col]])
    )
  

  map_data <- shapefile %>%
    left_join(rank_obj$table, by = admin_col)
  

  map_colors <- c(high_color, middle_color, low_color)
  names(map_colors) <- c(high_color, middle_color, low_color)
  
  map_labels <- c("High", "Middle", "Low")
  names(map_labels) <- c(high_color, middle_color, low_color)
  

  p_map <- ggplot(data = map_data) +
    geom_sf(aes(fill = color_group), color = "black", lwd = 0.1) +
    scale_fill_manual(
      values = map_colors,
      labels = map_labels,
      name   = "",
      na.value = "grey80"
    ) +
    labs(title = title) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      panel.grid   = element_blank()
    )
  
  return(p_map)
}


# **********************
# 7. map plots-------
# **********************

map_smth2 <- rank_map(
  rank_obj    = rank_smth2,
  shapefile   = poly.adm2,
  admin_col   = "admin2.name.full",
  high_color   = high_color,
  middle_color = middle_color,
  low_color    = low_color,
  title        = ""
)

map_smth2_fix <- rank_map(
  rank_obj    = rank_smth2_fix,
  shapefile   = poly.adm2,
  admin_col   = "admin2.name.full",
  high_color   = high_color,
  middle_color = middle_color,
  low_color    = low_color,
  title        = ""
)


# **********************
# 8. save plots-------
# **********************


# 1) Fay-Herriot ranking plots
ggplot2::ggsave(
  filename = here::here("figures", "dhs", "rank_smth2.png"),
  plot = rank_smth2$plot,
  width = 11,
  height = 16,
  dpi = 300
)


# 2) Fixed Fay-Herriot ranking plot
ggplot2::ggsave(
  filename = here::here("figures", "dhs", "rank_smth2_fix.png"),
  plot = rank_smth2_fix$plot,
  width = 11,
  height = 16,
  dpi = 300
)


# 3) FH model map plot
ggplot2::ggsave(
  filename = here::here("figures", "dhs", "map_smth2.png"),
  plot = map_smth2,
  width = 8,
  height = 8,
  dpi = 300
)


# 4) FH model fix map plot
ggplot2::ggsave(
  filename = here::here("figures", "dhs", "map_smth2_fix.png"),
  plot = map_smth2_fix,
  width = 8,
  height = 8,
  dpi = 300
)



