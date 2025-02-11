################################################################################
# Set up environment

rm(list = ls(all = TRUE))

# Load packages
source("R/utils/load_packages.R")

# Fix default ggsave background
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
theme_set(theme_minimal())

# Set region and time period of interest
cc = list(name = "TEST", proj = "EPSG:4263",
       start = "2019-01-01", end = "2021-01-01",
       type = "WPV1", cluster = NULL, eg = NULL)
# cc = list(name = "NGA", proj = "EPSG:4263",
#        start = "2014-09-01",end = "2016-08-01",
#        type = "WPV1", cluster = NULL, eg = NULL) # clust = "R4B5C5B2B3"
# cc = list(name = "NGA", proj = "EPSG:4263",
#        start = "2016-10-01", end = "2020-09-01",
#        type = "WPV1", cluster = NULL, eg = NULL)

################################################################################
################################################################################
