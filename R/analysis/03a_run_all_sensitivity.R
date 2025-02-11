################################################################################
################################################################################
# Time-rolling estimation of local risk and surveillance sensitivity
# *Sensitivity analyses*
# 
# + Primary analysis repeated across a given range of values for 
#   - ES catchment size
#   - Design prevalence
################################################################################
# Set up environment, directories 

rm(list = ls(all = TRUE))

# Load packages
source("R/utils/load_packages.R")

cc = c("TEST","EPSG:4263","2019-01-01","2021-01-01",NULL)
# cc = c("NGA","EPSG:4263","2014-09-01","2016-08-01",NULL)
cc_name = tolower(cc[1])

# Local projection
proj_local <- cc[2]

# Specify virus type to analyse
type = "WPV1"

# Specify genetic cluster to analyse
clust = cc[5]

# Inputs
indir <- paste0("inputs/",cc_name)

# Outputs
outdir <- ifelse(!is.na(clust),
                 here("results",cc_name,clust,substr(cc[3],1,4),"sensitivity analyses"),
                 here("results",cc_name,substr(cc[3],1,4),"sensitivity analyses"))

if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

# Fix default ggsave background
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

theme_set(theme_minimal(base_size = 16))

################################################################################
# Load data

# District data
data <- readRDS(here(indir, "data.rds"))

# Radiation matrix
R <- readRDS(here(indir, "R.rds"))

# AFP cases
afp <- readRDS(here(indir, "afp.rds")) 

# ES 
es <- readRDS(here(indir, "es.rds"))

# District shapes
shape2 <- readRDS(here(indir, "shape2.rds"))

# Population raster
poprast <- rast(here(indir,paste0("worldpop_raster_proj_",cc_name,".tif"))) %>% 
  project(poprast, proj_local)

# Detection sensitivities
AFPvals <- readxl::read_excel(here("inputs","AFPSurvValues.xlsx"))
ESvals <- readxl::read_excel(here("inputs","ESSurvValues.xlsx"))

# Source function
source("R/utils/run_all_fcn.R")
source("R/utils/estimate_elimination_fcn.R")

date_to_period = function(x){
  year(x) + (month(x)-1)/12
}

period_to_date = function(x,format='%Y-%b'){
  yr = floor(x)
  mo = round(12*(x-yr),0)+1
  d = '1'
  format(ymd(paste(yr,mo,d,sep='-')),format)
}

scale_val <- function(x){
  (x - min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))
}

################################################################################

# Number of replicates
Iter = 1000

## - Parameters to vary - ##
# Switch out commented code in order to vary catchment size or design prevalence

## Design prevalence (district level monthly rate)
# Fixed:
Dprev = 1/(100000/12)
# Varying:
# Dprev_range = 1/(c(1000,10000,100000)/12)
# grp_var = "dprev"
# grp_nm = "Design prevalence"
# grp_labs <- c("1/1,000","1/10,000","1/100,000")

## Catchment size (radius around point location of sampling site)
# Fixed:
# buffer_km = 5
# Varying:
buffer_range = c(2,5,10)
grp_var = "buff"
grp_nm = "Catchment radius"
grp_labs <- c("2km","5km","10km")

## To vary design prevalence:
# out <- lapply(Dprev_range,
#               function(Dprev){

## To vary catchment size:
out <- lapply(buffer_range,
              function(buffer_km){

                # Calculate rolling risk and sensitivity per month during period of interest
                tp_m <- c(seq(ymd(cc[3]),ymd(cc[4]), by = "month"))
                lapply(tp_m,
                       function(tdy){
                         run_FFI(tdy,
                                 data = data,
                                 afp = afp,
                                 es = es,
                                 R = R,
                                 shapes = shape2,
                                 poprast = poprast,
                                 proj_local = proj_local,
                                 buffer_km = buffer_km,
                                 return_catchment_shapes = F,
                                 Dprev = Dprev,
                                 Iter = Iter,
                                 equal_risk = F,
                                 est_FFI = F,
                                 M = 2)
                       }) -> out_mth
                
                out_all <- setNames(out_mth, tp_m)
                
               })

# ---------------------------------------------------------------------------- #
# SENSITIVITY ANALYSES

# Identify districts with any coverage
tp <- sapply(out[[1]], function(x) as.character(x$sens$sens_period[1]))
tp_f <- format(ymd(tp), "%b-%Y") 

data_df_setup <- function(obj){
  lapply(obj, 
         function(x) x$sens$data) %>% 
    setNames(tp) %>% 
    bind_rows(.id = "tdy") %>% 
    mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                          levels = tp_f)) -> out
  return(out)
}
data_all <- lapply(out, data_df_setup) %>%
  setNames(grp_labs) %>%
  bind_rows(.id = "grp") %>%
  mutate(grp = factor(grp, levels = grp_labs))

dist_with_es <- data_all %>% 
  arrange(tdy,tdy_f) %>%
  select(grp, es_coverage) %>% 
  group_by(grp) %>% group_split() %>% 
  lapply(function(x) pull(x,es_coverage))

# Extract sensitivity and pneg estimates and combine
sens_df_setup <- function(obj, pneg = F){
  
  lapply(obj, 
         function(x) x$sens$sens_df) %>% 
    setNames(tp) %>% 
    bind_rows(.id = "tdy") %>% 
    mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                          levels = tp_f)) -> out
  
  if(pneg == T){
    # Extract Pneg draws for ENV to recalculate sensitivity in districts with coverage
    lapply(obj,
           function(x) as_tibble(t(x$sens$sens_dist_draws$Pneg_ENVi))) %>%
      setNames(tp) %>%
      bind_rows(.id = "tdy") %>%
      mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                            levels = tp_f)) %>% 
      select(tdy, tdy_f, everything()) -> pneg_ENV
    return(pneg_ENV)
 
    }else{return(out)}
}

sens_all <- lapply(out, sens_df_setup) %>%
  setNames(grp_labs) %>%
  bind_rows(.id = "grp") %>%
  mutate(grp = factor(grp, levels = grp_labs))

saveRDS(sens_all, here(outdir, paste0("sens_by",grp_var,".rds"))) 

# Extract pneg estimates, filter to only LGAs with >X% ES coverage, and recalculate sensitivity overall
pneg_env_all <-  lapply(out, sens_df_setup, pneg = T) %>%
  setNames(grp_labs)

calc_sens_subset <- function(pneg, dists){
  dists %>% 
    bind_cols(pneg) %>%
    filter(es_coverage > 0.25) %>%
    pivot_longer(-1:-4, values_to = "Pneg_ENVi", names_to = "Iter") %>% #ggplot(aes(x = Pneg_ENVi, fill = (es_coverage > 0.25))) + geom_histogram(col = NA, position = "dodge",alpha = 0.5) + scale_y_continuous(trans = "sqrt") 
    group_by(tdy, tdy_f,Iter) %>% 
    summarise(SSe_ENV = 1-prod(Pneg_ENVi)) %>% 
    ungroup() %>% 
    return()
}

sens_ENV_subset <- purrr::map2(pneg_env_all, dist_with_es, calc_sens_subset) %>%
  setNames(grp_labs) %>%
  bind_rows(.id = "grp") %>%
  mutate(grp = factor(grp, levels = grp_labs),
         cat = ">25% coverage")

# ---------------------------------------------------------------------------- #
# Plot

# Catchment radius
sens_all %>%
  # select(grp, tdy, tdy_f, SSe_ENV) %>% 
  # mutate(cat = "All") %>% 
  # bind_rows(sens_ENV_subset) %>%
  ggplot(aes(grp, SSe_ENV, col = grp)) + #, lty = cat
  geom_boxplot() +
  guides(col = "none") +
  scale_y_continuous(trans = "log10") +
  scale_colour_viridis_d(option = "rocket", begin = 0.3, end = 0.8, direction = 1) +
  # scale_colour_manual(values = pal[-2]) +
  theme(legend.position = "bottom") +
  labs(x = grp_nm,y = "Sensitivity",
       col = grp_nm,
       lty = NULL
       # title = "Estimated monthly surveillance sensitivity",
       # subtitle = "By assumed catchment radius around active ES sites"
  ) -> p_bybuff
p_bybuff
ggsave(here(outdir,"figures/sensitivity analyses","sens_bycomp_bybuff.png"), p_bybuff, height=5,width=6)
# ggsave(here(outdir,"figures/sensitivity analyses","sens_bycomp_bybuff_bycovg.png"), p_bybuff, height=5,width=6)

# National ES coverage by assumed catchment radius
data_all %>% 
  # Sum across LGAs for each tp
  group_by(grp, tdy) %>% 
  summarise(catchment_pop = sum(catchment_pop),
            total_pop = sum(total_pop), 
            es_coverage = catchment_pop/total_pop) %>% 
  ungroup() -> nat_covg 

nat_covg %>% 
  ggplot(aes(grp, es_coverage)) +
  geom_boxplot() +
  # scale_y_continuous(trans = "log10") +
  labs(x = "Assumed catchment radius", y = "National coverage")
ggsave(here(outdir,"figures/sensitivity analyses","es_covg_bybuff.png"), height=4,width=5)

# Average across time
nat_covg %>% 
  group_by(grp) %>% 
  summarise(es_coverage = mean(es_coverage)) 
# 1 2km        0.0117
# 2 5km        0.0291
# 3 10km       0.0425

# Sensitivity across all LGAs with/without threshold coverage
sens_all %>%
  select(grp, tdy, tdy_f, SSe_ENV) %>% 
  mutate(cat = "All") %>% 
  bind_rows(sens_ENV_subset) -> sens_bycovg
  
sens_bycovggrp %>%  
  group_by(grp, cat) %>% 
  summarise(mean_sens = mean(SSe_ENV),
            iqr_sens = paste(round(quantile(SSe_ENV, c(0.25,0.75)),4), collapse = ","))

#   grp   cat           mean_sens iqr_sens     
# 1 2km   >25% coverage    0.0122 0.0112,0.0131
# 2 2km   All              0.0210 0.0196,0.0222

# 3 5km   >25% coverage    0.0213 0.0205,0.0222
# 4 5km   All              0.0280 0.0264,0.0297

# 5 10km  >25% coverage    0.0329 0.0301,0.0357
# 6 10km  All              0.0343 0.0317,0.0371

# ---------------------------------------------------------------------------- #

# Design prevalence
sens_all %>%
  pivot_longer(c("SSe_AFP","SSe_ENV")) %>%
  mutate(name = factor(name, labels = c("AFP","ENV"))) %>%
  ggplot(aes(grp, value, col = grp)) + 
  geom_boxplot() +
  facet_wrap(~name, ncol = 2, scales = "free_y") +
  guides(col = "none") +
  # scale_y_continuous(trans = "log10") +
  scale_colour_viridis_d(option = "rocket", begin = 0.3, end = 0.8, direction = 1) +
  # scale_colour_manual(values = pal[-2]) +
  labs(x = NULL,#"Design prevalence",
       y = "Sensitivity",
       col = "Design prevalence"
       # title = "Estimated monthly surveillance sensitivity",
       # subtitle = "By assumed design prevalence of infection (annual rate) among the total population"
  ) -> p_bydprev
p_bydprev
ggsave(here(outdir,"figures/sensitivity analyses","sens_bycomp_bydprev.png"), p_bydprev, height=5,width=10)

################################################################################
# FFI

prior = c(30,30)
sens_all %>% 
  group_by(grp) %>% 
  group_split() %>% 
  lapply(function(sens_all){elim_est_tv(sens_all, 
                                        intro = T,
                                        prior_parms = prior)}) -> ffi_by_var

lapply(ffi_by_var, 
       function(x) x$free_AE_sims) %>% 
  setNames(grp_labs) %>% 
  bind_rows(.id = "grp") %>% 
  mutate(t = period_to_date(time),
         grp = factor(grp, levels = grp_labs)) -> ffi_compare

saveRDS(ffi_compare, here::here(outdir, paste0("ffi_compare_",grp_var,".rds")))

# Identify time point (if any) at which confidence criterion for elimination is reached
elim_criterion = 0.95
ffi_compare %>% 
  mutate(max_time = max(time)) %>% 
  group_by(grp,max_time) %>% 
  filter(med > elim_criterion) %>%
  slice_min(time) %>% 
  ungroup() -> ffi_thresholds

# Plot
pd = position_dodge2(width = 0.05)
ffi_compare %>% 
  ggplot(aes(x=time,y=med)) + 
  # Formatting
  geom_hline(aes(yintercept=0.95),linetype="dashed",col="grey50") +
  geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
  # geom_hline(aes(yintercept=1),linetype="solid",col="grey10") +
  geom_text(aes(x = min(time),y = 0.93), label = "0.95",col="grey50", cex = 2) +
  geom_text(aes(x = min(time),y = 0.97), label = "0.99",col="red", cex = 2) +
  # Elimination thresholds
  geom_rect(data = ffi_thresholds,
            aes(xmin = time,
                xmax = round(max_time,1),
                ymin = 0,
                ymax = 1,
                fill = grp),
            alpha = 0.1) +
  # Monthly FFI probabilities
  geom_pointrange(aes(ymin=lwr, ymax=upr,
                      col = grp),
                  position = pd) +
  scale_colour_viridis_d(option = "rocket", begin = 0.3, end = 0.8) +
  scale_fill_viridis_d(option = "rocket", begin = 0.3, end = 0.8) +
  scale_x_continuous(breaks = unique(round(ffi_compare$time))) +
  theme(legend.position = "bottom") +
  guides(fill = "none") +
  labs(x = "Month", 
       y = "Pr(Infection free)",
       col = grp_nm) + 
  ylim(0,1) -> p_ffi_compare

p_ffi_compare
ggsave(here(outdir,
            "figures/sensitivity analyses",
            paste0("compare_ffi_",grp_var,"_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
       p_ffi_compare, height=6,width=10)

# p_comb <- p_bybuff + p_ffi_compare + plot_layout(widths = c(1,2)) + plot_annotation(tag_levels = "A")
p_comb <- p_bydprev / p_ffi_compare + plot_layout(heights = c(2,3)) +plot_annotation(tag_levels = "A")
p_comb

ggsave(here(outdir,
            "figures/sensitivity analyses",
            paste0("comb_compare_ffi_",grp_var,"_",year(min(sens_all$tdy)),"-",year(max(sens_all$tdy)),".png")),
       p_comb,
       # height=5, width = 10
       height=9,width=10
       )

################################################################################
################################################################################