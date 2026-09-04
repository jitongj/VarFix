# *******************
# date: 2026.01.21
# task: This script is for downloading pop data from worldpop and generating k0-5
# author: Jitong Jiang
# ********************

rm(list = ls())
# ENTER COUNTRY OF INTEREST 
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Zambia'
gadm.abbrev <- 'zmb'
survey_year <- '2018'   ### if the survey spans two years, use the first year
survey <- paste0(country,'_',survey_year)
pop.abbrev <- tolower(gadm.abbrev)



# Load libraries and info 
library(dplyr)
library(sf)
library(terra)
library(rdhs)
library(caret)

subpop_dir <- here::here("data", "subpop")
dir.create(subpop_dir, recursive = TRUE, showWarnings = FALSE)


ages <- c(0, 1)

for (age in ages) {
  
  female_pop_file <- file.path(
    subpop_dir,
    paste0(pop.abbrev, "_f_", age, "_", survey_year, "_1km.tif")
  )
  
  if (!file.exists(female_pop_file)) {
    
    year_url <- if (survey_year < 2000) 2000 else survey_year
    
    url <- paste0(
      "https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_1km/unconstrained/",
      year_url, "/", toupper(pop.abbrev), "/",
      pop.abbrev, "_f_", age, "_", year_url, "_1km.tif"
    )
    
    download.file(url, female_pop_file, method = "libcurl", mode = "wb")
  }
  
  male_pop_file <- file.path(
    subpop_dir,
    paste0(pop.abbrev, "_m_", age, "_", survey_year, "_1km.tif")
  )
  
  if (!file.exists(male_pop_file)) {
    
    year_url <- if (survey_year < 2000) 2000 else survey_year
    
    url <- paste0(
      "https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_1km/unconstrained/",
      year_url, "/", toupper(pop.abbrev), "/",
      pop.abbrev, "_m_", age, "_", year_url, "_1km.tif"
    )
    
    download.file(url, male_pop_file, method = "libcurl", mode = "wb")
  }
}


kid_0_5_filename <- file.path(
  subpop_dir,
  paste0(pop.abbrev, "_k0_5_", survey_year, "_1km.tif")
)

if (!file.exists(kid_0_5_filename)) {
  
  ages <- c(0, 1)
  kid_0_5_raster <- NULL
  
  for (age in ages) {
    
    female_raster <- terra::rast(
      file.path(subpop_dir, paste0(pop.abbrev, "_f_", age, "_", survey_year, "_1km.tif"))
    )
    male_raster <- terra::rast(
      file.path(subpop_dir, paste0(pop.abbrev, "_m_", age, "_", survey_year, "_1km.tif"))
    )
    
    kid_raster <- female_raster + male_raster
    
    kid_0_5_raster <- if (is.null(kid_0_5_raster)) {
      kid_raster
    } else {
      kid_0_5_raster + kid_raster
    }
  }
  
  terra::writeRaster(kid_0_5_raster, kid_0_5_filename)
}
