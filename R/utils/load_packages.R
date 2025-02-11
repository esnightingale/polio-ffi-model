################################################################################
# Install/load all required packages
################################################################################

if(!require("pacman", character.only = TRUE)){
  install.packages("pacman")
}

pacman::p_load(here,
                 tidyverse,
                 sf,
                 lubridate,
                 scales,
                 patchwork,
                 terra,
                 ggrepel,
                 ggtext,
                 prevalence,
                 units,
                 zoo,
               install = T)

# Helper functions
date_to_period = function(x){
  year(x) + (month(x)-1)/12
}

period_to_date = function(x,format='%Y-%b'){
  yr = floor(x)
  mo = round(12*(x-yr),0)+1
  d = '1'
  format(ymd(paste(yr,mo,d,sep='-')),format)
}

# Fix default ggsave background
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

theme_set(theme_minimal(base_size = 14)) 

################################################################################
################################################################################
