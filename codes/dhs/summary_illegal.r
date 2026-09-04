# *******************
# date: 2026.01.21
# task: This script is for table of illegal admin2 areas
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


data_urban <- data[data$strata == "urban", ]
data_rural <- data[data$strata == "rural", ]

res_ad0 <- directEST(data = data,
                     cluster.info = cluster.info,
                     admin = 0)
res_ad0_urban <- directEST(data = data_urban,
                           cluster.info = cluster.info,
                           admin = 0)
res_ad0_rural <- directEST(data = data_rural,
                           cluster.info = cluster.info,
                           admin = 0)
# **********************
# 3. national weight-------
# **********************

nat_weight <- tibble(
  strata = c("urban", "rural"),
  nat_avg_total_weight = c(16662007, 18063987)
)

# 2. combine with direct.est
national_df <- nat_weight %>%
  mutate(
    direct.est = c(
      res_ad0_urban$res.natl$direct.est,
      res_ad0_rural$res.natl$direct.est
    )
  )

national_df



# **********************
# 4. illegal data -------
# **********************

illegal_index <- tibble::tibble(
  Index = 1:24,
  admin2 = c(
    "Chitambo", "Kapiri Mposhi", "Luano", "Ngabwe",
    "Chingola",
    "Chadiza", "Chasefu", "Mambwe", "Vubwi",
    "Chembe", "Chipili", "Milengi",
    "Lavushimanda",
    "Chavuma", "Ikelenge", "Mufumbwe", "Mushindano",
    "Kaputa", "Mporokoso",
    "Kaoma", "Mitete", "Mwandi", "Nalolo", "Sioma"
  ),
  admin1 = c(
    rep("Central", 4),
    "Copperbelt",
    rep("Eastern", 4),
    rep("Luapula", 3),
    "Muchinga",
    rep("North-Western", 4),
    rep("Northern", 2),
    rep("Western", 5)
  )
) %>%
  mutate(
    admin2.name.full = paste0(admin1, "_", admin2)
  )




# **********************
# 3. FH estimate results-------
# **********************

bad_admin2 <- res_ad2_fix$fixed_areas

country_shp_analysis <- readRDS(here::here("data/Zambia/country_shp_analysis.rds"))

data0 <- data %>% filter(!is.na(value))

myData_tmp <- data0 %>%
  group_by(cluster) %>%   # Group by cluster
  dplyr::summarise(       # Compute within-cluster summary statistics
    strata = unique(strata),  # Unique strata value per cluster (assuming one strata per cluster)
    Ntrials = n(),            # Number of individuals in the cluster
    weight = unique(weight),
    value = sum(value),       # Total number of successes (sum of individual values)
    households_number = n_distinct(householdID) # Number of distinct households in the cluster
  )

# combine individual level data with cluster info
myData <-merge(myData_tmp, cluster.info[["data"]], by="cluster")

sampled_admin2_cluster <- myData %>%
  group_by(admin1.name, admin2.name, admin2.name.full) %>%
  summarise(sampled_admin2_cluster = n(), .groups = "drop")


# bad_admin2 is the list of admin2 areas returned by res_ad2_fix$fixed_areas
bad_admin2_list <- bad_admin2

# Extract raw cluster-level observations corresponding to bad_admin2 areas
bad_admin2_raw <- myData %>%
  filter(admin2.name.full %in% bad_admin2_list) %>%
  transmute(
    Admin1 = admin1.name,
    admin2.name.full = admin2.name.full,
    strata = strata,
    clusterID = cluster,
    outcome = value,
    Ntrials = Ntrials,
    weight = weight
  ) %>%
  arrange(admin2.name.full) 

bad_admin2_raw

# **********************
# 4. bad data points -------
# **********************

bad_admin2_with_nat <- bad_admin2_raw %>%
  left_join(national_df, by = "strata")

bad_admin2_final <- bad_admin2_with_nat %>%
  left_join(
    illegal_index %>% select(Index, admin2, admin2.name.full),
    by = "admin2.name.full"
  )

View(bad_admin2_final)

illegal_data <- bad_admin2_final %>%
  transmute(
    index           = Index,
    Admin2          = admin2,
    Admin1          = Admin1,
    `Urban/Rural`   = strata,      
    clusterID       = clusterID,
    n_trials        = Ntrials,
    outcome         = outcome,
    phantom_mean    = round(direct.est, 3),  # national direct.est
    phantom_weights = nat_avg_total_weight
  ) %>%
  arrange(index, Admin1, Admin2, `Urban/Rural`, clusterID)

illegal_data



# 
# write.csv(illegal_data,
#           file = "illegal_admin2_raw.csv",
#           row.names = FALSE)

