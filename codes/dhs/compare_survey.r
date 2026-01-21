# *******************
# date: 2026.01.21
# task: This script is for comparing unfix sas method with survey package, 
# author: Jitong Jiang
# ********************

source(here::here("codes", "directEST_1030_national.R"))

library(dplyr)
library(purrr)
library(ggrepel)
library(kableExtra)
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
library(png)
library(grid) 
library(sf)
library(viridis)
library(gridExtra)


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
# 2. survey design   -------
# **********************

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


data_urban <- data[data$strata == "urban", ]
data_rural <- data[data$strata == "rural", ]

data$v022 = paste(data$v024, "-", data$strata)
data$v023 = paste(data$v024, "-", data$strata)

modt <- left_join(data, cluster.info$data, by = "cluster")
modt <- modt[!(is.na(modt$admin2.name)), ]
modt <- modt[!(is.na(modt$value)), ]
modt$admin1_strata <- paste0(modt$admin1.name,"_",modt$strata)
# calculate L_ih
modt <- modt %>%
  group_by(v023) %>%
  mutate(L_ih = n_distinct(cluster)) %>%
  ungroup()


modt_new = modt
modt_new$admin1_strata = paste0(modt_new$admin1.name,"_",modt_new$strata)
modt_new$admin1_strata <- as.factor(modt_new$admin1_strata)

modt_new$admin2_strata = paste0(modt_new$admin2.name.full,"_",modt_new$strata)
modt_new$admin2_strata <- as.factor(modt_new$admin2_strata)

DHSdesign<-survey::svydesign(id=modt_new$cluster, strata=modt_new$v023, weights=modt_new$weight, data=modt_new)


# **********************
## 2.1. stratified survey direct estimate -------
# **********************

# urban/rural estimate

### admin1
mean_strat_admin1 <- survey::svyby(~value, ~admin1_strata, DHSdesign, svymean, vartype=c("se"), 
                                   drop.empty.groups=FALSE)
mean_strat_admin1$strata <- ifelse(grepl("urban", mean_strat_admin1$admin1_strata), "urban",
                                   ifelse(grepl("rural", mean_strat_admin1$admin1_strata), "rural", NA))
mean_strat_admin1$admin1.name <- sub("_(urban|rural)$", "", mean_strat_admin1$admin1_strata)

mean_strat_admin1$direct.var = mean_strat_admin1$se^2

direct_admin1_urban <- mean_strat_admin1 %>%
  filter(strata == "urban") %>%
  rename(direct.est = value)

direct_admin1_rural <- mean_strat_admin1 %>%
  filter(strata == "rural") %>%
  rename(direct.est = value)


### admin2
mean_strat_admin2 <- survey::svyby(~value, ~admin2_strata, DHSdesign, svymean, vartype=c("se"), 
                                   drop.empty.groups=FALSE)
mean_strat_admin2$strata <- ifelse(grepl("urban", mean_strat_admin2$admin2_strata), "urban",
                                   ifelse(grepl("rural", mean_strat_admin2$admin2_strata), "rural", NA))
mean_strat_admin2$admin2.name.full <- sub("_(urban|rural)$", "", mean_strat_admin2$admin2_strata)

mean_strat_admin2$direct.var = mean_strat_admin2$se^2

direct_admin2_urban <- mean_strat_admin2 %>%
  filter(strata == "urban") %>%
  rename(direct.est = value) %>%
  mutate(
    admin2.name.full = as.character(admin2.name.full),
    admin1.name = sub("_.*", "", admin2.name.full)
  )

direct_admin2_rural <- mean_strat_admin2 %>%
  filter(strata == "rural") %>%
  rename(direct.est = value) %>%
  mutate(
    admin2.name.full = as.character(admin2.name.full),
    admin1.name = sub("_.*", "", admin2.name.full)
  )


# **********************
## 2.2. unstratified survey direct estimate -------
# **********************

### admin1
mean_admin1 <- survey::svyby(~value, ~admin1.name, DHSdesign, svymean, vartype=c("se"), 
                             drop.empty.groups=FALSE)
mean_admin1$direct.var = mean_admin1$se^2
direct_admin1 = mean_admin1

direct_admin1 <- direct_admin1%>%
  as.data.frame() %>%
  filter(is.na(value) == FALSE) %>%
  distinct()

direct_admin1 <- direct_admin1 %>%
  as.data.frame() %>%
  filter(!is.na(value)) %>%
  rename(direct.est = value) %>%
  distinct()

### admin2
mean_admin2 <- survey::svyby(~value, ~admin2.name.full, DHSdesign, svymean, vartype=c("se"), 
                             drop.empty.groups=FALSE)
mean_admin2$direct.var = mean_admin2$se^2
direct_admin2 = mean_admin2

direct_admin2 <- direct_admin2 %>%
  as.data.frame() %>%
  filter(!is.na(value)) %>%
  rename(direct.est = value) %>%
  distinct() %>%
  mutate(
    admin2.name.full = as.character(admin2.name.full), 
    admin1.name = sub("_.*", "", admin2.name.full)  
  )


# **********************
# 3. unfix sas results   -------
# **********************

# DHS data
res_ad1_sas <- directEST_1030(data = data,
                             cluster.info = cluster.info,
                             admin = 1)

# 1.2 compare at admin2 level
res_ad2_sas <- directEST_1030(data = data,
                             cluster.info = cluster.info,
                             admin = 2, var.fix=FALSE)


# **********************
# 4. comapre sas and survey   -------
# **********************

cmp_adm1_est <- data.frame(
  survey_est = direct_admin1$direct.est,
  sas_est    = res_ad1_sas$res.admin1$direct.est,
  admin1     = direct_admin1$admin1.name
)

cmp_adm1_var <- data.frame(
  survey_var = direct_admin1$direct.var,
  sas_var    = res_ad1_sas$res.admin1$direct.var,
  admin1     = direct_admin1$admin1.name
)

cmp_adm2_est <- data.frame(
  survey_est = direct_admin2$direct.est,
  sas_est    = res_ad2_sas$res.admin2$direct.est,
  admin2     = direct_admin2$admin2.name.full
)

cmp_adm2_var <- data.frame(
  survey_var = direct_admin2$direct.var,
  sas_var    = res_ad2_sas$res.admin2$direct.var,
  admin2     = direct_admin2$admin2.name.full
)


# **********************
# 5. plots   -------
# **********************
p1 <- ggplot(cmp_adm1_est, aes(x = survey_est, y = sas_est)) +
  geom_point(size = 2, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    x = "survey",
    y = "SAS",
    title = "Admin1: direct.est"
  ) +
  theme_minimal() +
  coord_equal()

p2 <- ggplot(cmp_adm1_var, aes(x = survey_var, y = sas_var)) +
  geom_point(size = 2, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    x = "survey",
    y = "SAS",
    title = "Admin1: direct.var"
  ) +
  theme_minimal() +
  coord_equal()

p3 <- ggplot(cmp_adm2_est, aes(x = survey_est, y = sas_est)) +
  geom_point(size = 2, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    x = "survey",
    y = "SAS",
    title = "Admin2: direct.est"
  ) +
  theme_minimal() +
  coord_equal()

p4 <- ggplot(cmp_adm2_var, aes(x = survey_var, y = sas_var)) +
  geom_point(size = 2, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    x = "survey",
    y = "SAS",
    title = "Admin2: direct.var"
  ) +
  theme_minimal() +
  coord_equal()


sas_vs_survey_unfix <- (p1 + p2) / (p3 + p4)

saveRDS(
  sas_vs_survey_unfix,
  file = here::here("figures","dhs" , "survey_compare", "sas_vs_survey_unfix.rds")
)

