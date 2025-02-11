################################################################################
################################################################################
# Analysis data setup
# + Subsetting to region of interest
# + Selecting relevant variations
# + Defining new variables
#    - WPV +/-, population, immunity, space/time indices
# + Defining matching district IDs
# + Defining radiation matrix across districts
#
# -> Save outputs for analysis
#
# Note: This script does not need to be run if using pre-defined test inputs
################################################################################

library(tidyverse)
library(sf)
library(zoo)
library(janitor)
library(units)
library(Matrix)
library(config)
library(here)
library(lubridate)
library(terra)

rm(list = ls(all = TRUE))  # clear all data

source("R/utils/setup_env.R")
cc_name <- tolower(paste0(cc$name, collapse = "_"))

# Local projection
proj_local <- cc$proj

################################################################################
# Set up environment, directories

set_here()

# Inputs
filename_tsir_data = here("inputs","tsir_data.csv")

indir <- "inputs"
filename_afp = here(indir,"linelist_afp_clean.csv")

filename_es = here(indir,"linelist_es_clean.csv")

filename_npafp = here(indir,"npafp_clean.csv")

# Outputs
outdir <- paste0(indir,cc_name)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

# ---------------------------------------------------------------------------- #
# Additional helper functions

period_to_date = function(x,format='%Y-%b'){
  yr = floor(x)
  mo = round(12*(x-yr),0)+1
  d = '1'
  format(ymd(paste(yr,mo,d,sep='-')),format)
}

date_to_period = function(x){
  year(x) + (month(x)-1)/12
}

################################################################################
# Read in Data #

# Overall observations per district, period and serotype
# - immunity, pop size, routine imm (dpt3) covg, cases, ES samples/positives
data = read_csv(filename_tsir_data)

# AFP case records
afp = read_csv(filename_afp)

# Environmental surveillance sample records
es = read_csv(filename_es)

# NPAFP for denominators
npafp <- read_csv(filename_npafp)

# District shapefiles
shape2 <- readRDS(here("inputs","shape2.rds"))
shape1 <- readRDS(here("inputs","shape1.rds"))
shape0 <- readRDS(here("inputs","shape0.rds"))

# Population raster
poprast <- rast(here("inputs",cc_name,paste0(cc_name,"_ppp_2020_UNadj.tif")))

################################################################################
# Set up data for analysis

# Add iso3 code to other datasets
afp = afp %>%
  left_join(shape0 %>%
              st_drop_geometry() %>%
              distinct(iso_3_code,adm0_name))

es = es %>%
  dplyr::select(iso_3_code:yr,
                site_id, site_name, site_name:guid,
                collection_date, travel_time, sample_condition, final_class,
                npev, virus_type_s, virus_cluster_s, emergence_group_s,
                coord_imp, x, y)

################################################################################
# Other variables needed for model

# Indices
period_index = tibble(period = sort(unique(data$period))) %>%
  mutate(i = 1:n())
area_index = tibble(guid = unique(data$guid)) %>%
  mutate(j=1:n())

data = left_join(data,period_index) %>%
  left_join(area_index)
data = data %>%
  arrange(serotype,j,i) %>%
  # susceptibility per district/month
  mutate(s = pmin(1,1-immunity))

# Add U5 populations from district data
shape2 = left_join(shape2,data %>%
                     distinct(guid,population_u5,j))
shape2 = shape2 %>% arrange(j)

# ---------------------------------------------------------------------------- #
# Subset to country/countries of interest

data_sub = data %>% filter(iso_3_code %in% cc)

# Also exclude AFP among over-15s
afp_sub = afp %>% filter(iso_3_code %in% cc) %>% filter(age_m < 15*12)

es_sub = es %>% filter(iso_3_code %in% cc)

npafp_sub = npafp %>% filter(iso_3_code %in% cc)

t <- as.Date("2023-06-01")
shape2_sub = shape2 %>% filter(enddate > t,
                           iso_3_code %in% cc)
shape1_sub = shape1 %>% filter(enddate > t,
                           iso_3_code %in% cc)
shape0_sub = shape0 %>% filter(enddate > t,
                           iso_3_code %in% cc)

# ---------------------------------------------------------------------------- #
# Population denominators for Nigeria ONLY

# adm2 population under 15yo
# + sourced from HDX / Nigeria data grid: https://data.humdata.org/dataset/cod-ps-nga?force_layout=desktop

if(cc == "NGA"){
  pop_hdx <- read_csv("inputs/nga/nga_admpop_adm2_2020_corrected.csv") %>%
    mutate(pop_u5 = T_00_04,
           pop_u15 = T_00_04 + T_05_09 + T_10_14,
           pop_hdx = T_TL) %>%
    select(ADM0_NAME,ADM1_NAME_MATCHED,ADM2_NAME_MATCHED, pop_hdx, pop_u5, pop_u15) %>%
    rename_with(tolower) %>%
    mutate(across(adm0_name:adm2_name_matched, toupper),
           adm2_name_matched = gsub("  ", " ",
                                    gsub("-"," ",adm2_name_matched)))

  shape2 <- shape2 %>%
    mutate(adm2_name = gsub("  ", " ",
                            gsub("-"," ",adm2_name))) %>%
    left_join(pop_hdx, by = c("adm0_name",
                              "adm1_name" = "adm1_name_matched",
                              "adm2_name" = "adm2_name_matched"))
}else{
  shape2 <- mutate(shape2, pop_hdx = NA, pop_u5 = NA, pop_u15 = NA)
}

# shape2_ %>%
#   filter(iso_3_code=="NGA" & is.na(pop_hdx)) %>%
#   select(adm1_name, adm2_name) %>%
#   distinct() -> nonmatch

# ---------------------------------------------------------------------------- #
# Locations and catchment of all ES sampling sites

poprast <- project(poprast, proj_local)

# Estimate total 5km catchment population for each site
# (this is recalculated in main monthly analysis for each actively sampling site,
# counting only catchment within the same district)
es_sub_catch5 <- es_sub %>%
  st_as_sf(coords = c("x","y"), crs = 4326) %>%
  st_transform(proj_local) %>%
  st_buffer(dist = units::set_units(5, "km"))

es_sub$catchment_pop_5k <- exactextractr::exact_extract(poprast,
                                                        es_sub_catch5,
                                                        fun = "sum")

# Extract unique sites
es_sub %>%
  select(iso_3_code:site_name, x, y, coord_imp, guid, catchment_pop_5k) %>%
  distinct() -> es_sites

# ---------------------------------------------------------------------------- #
# Reproject spatial data

shape2_sub <- st_transform(shape2_sub, proj_local)
shape1_sub <- st_transform(shape1_sub, proj_local)
shape0_sub <- st_transform(shape0_sub, proj_local)

# Extract total district populations from raster
shape2_sub$total_pop <- exactextractr::exact_extract(poprast,
                                                     shape2_sub,
                                                     fun = "sum")

data_sub <- left_join(data_sub, npafp_sub %>% select(guid, period, pop_u15_polisdenom))

summary(data_sub$pop_u15_polisdenom)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 3773   54379  101295  163528  154966 5311744    7528

# *Nigeria only*

if(cc == "nga"){
  # Check missing u15 pop -> impute based on u5 pop and ratio
  summary(shape2_sub$pop_u15)

  # One missing in NGA - impute with average ratio from u5 population
  r <- median(shape2_sub$pop_u15/shape2_sub$population_u5, na.rm = T) # ~2
  shape2_sub$pop_u15[is.na(shape2_sub$pop_u15)] <- shape2_sub$population_u5[is.na(shape2_sub$pop_u15)]*r
}

data_sub <- left_join(data_sub, shape2_sub %>% st_drop_geometry() %>% select(guid, total_pop, pop_hdx, pop_u5, pop_u15))

r <- median(data_sub$pop_u15_polisdenom/data_sub$total_pop, na.rm = T)
data_sub$pop_u15_polisdenom[is.na(data_sub$pop_u15_polisdenom)] <- data_sub$total_pop[is.na(data_sub$pop_u15_polisdenom)]*r

# ---------------------------------------------------------------------------- #
# Radiation matrix
# - Only needed for variable construction
# - Sparse matrix, same dim as shape2

R = radiation_matrix_construction(shape2_sub %>%
                                    st_drop_geometry %>%
                                    dplyr::select(center_lon,center_lat),
                                  shape2_sub$population_u5)
dim(R)

################################################################################
# Save analysis datasets

# District data
saveRDS(data_sub, here(outdir, "data.rds"))

# Radiation matrix
saveRDS(R, here(outdir, "R.rds"))

# AFP cases
saveRDS(afp_sub, here(outdir, "afp.rds"))

# ES
saveRDS(es_sub, here(outdir, "es.rds"))
saveRDS(es_sites, here(outdir, "es_sites.rds"))

# District shapefiles
saveRDS(shape2_sub, here(outdir, "shape2.rds"))
saveRDS(shape1_sub, here(outdir, "shape1.rds"))
saveRDS(shape0_sub, here(outdir, "shape0.rds"))

# Projected population raster
writeRaster(poprast,
            here(outdir,paste0("worldpop_raster_proj_",cc_name,".tif")),
                 overwrite=T)

################################################################################
################################################################################
