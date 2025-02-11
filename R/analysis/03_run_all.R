################################################################################
################################################################################
# Time-rolling estimation of local risk and surveillance sensitivity
#
# + Estimate sensitivity each month for a given set of admin units, according to
#   an assumed ES catchment size and design prevalence.
# + Outputs are summarised and saved in "results/.."
# + Surveillance effort and estimated sensitivity per month/admin unit are
#   also mapped and saved in "results/../figures"
################################################################################
# Set up environment, directories

# Set region and time period of interest
source("R/utils/setup_env.R")

cc_name <- tolower(paste0(cc$name, collapse = "_"))

# Local projection (for drawing buffers & estimating ES catchment)
proj_local = cc$proj

# Virus type to analyse
type = cc$type

# Specify genetic cluster to analyse
clust = cc$cluster

# Specify emergence group to analyse
eg = cc$eg

# Inputs
indir <- paste0("inputs/",cc_name)

# Outputs
if(!is.null(clust)){outdir <- here("results",cc_name,type,clust,substr(cc$start,1,4),)
}else if (!is.null(eg)){outdir <- here("results",cc_name,type,eg,substr(cc$start,1,4))
}else {outdir <- here("results",cc_name,type,substr(cc$start,1,4))}

if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

################################################################################
# Load data

# District data
data <- readRDS(here(indir, "data.rds")) %>%
  # Retain only for specified serotype
  filter(serotype == gsub("[^0-9.-]", "", type))

# Radiation matrix
R <- readRDS(here(indir, "R.rds"))

# AFP cases
afp <- readRDS(here(indir, "afp.rds"))

# ES
es <- readRDS(here(indir, "es.rds"))

# District shapes
shape2 <- readRDS(here(indir, "shape2.rds"))

# Population raster
poprast <- rast(here(indir,paste0("worldpop_raster_proj_",cc_name,".tif")))

################################################################################
# Set assumed detection sensitivities - depending on virus type

# For WPV1
if (type == "WPV1"){
  AFPvals <- data.frame(SurveillanceNode = c("AFPcase",       # Pr. infection develops into disease
                                             "AFPnotified",    # Pr. adequate stool sample collected
                                             "AFPtest"),      # Viral load in stool sample above LoD - will be affected by sample quality and timeliness of analysis
                        EstMean = c(190,0.9,0.97), #WPV1 values
                        EstLwr = c(250,0.6,0.95),
                        EstUpr = c(150,0.999,0.999))
}else if (type == "cVDPV2"){
  # For cVPDV2
  AFPvals <- data.frame(SurveillanceNode = c("AFPcase",
                                             "AFPnotified",
                                             "AFPtest"),
                        EstMean = c(2000,0.9,0.97),
                        EstLwr = c(2200, 0.6,0.95),
                        EstUpr = c(1800,0.999,0.999))
}

# For WPV1/cVPDV2
ESvals <- data.frame(SurveillanceNode = "EStest", # Virus load above LoD - will be affected by site factors, conditions of sample collection and timeliness of analysis
                      EstMean = 0.9,
                      EstLwr = 0.7,
                      EstUpr = 0.99)

################################################################################
# Source functions/plotting options

source("R/utils/run_all_fcn.R")
source("R/utils/mapping_fcns.R")
source("R/utils/helper_fcns.R")

figdir <- here(outdir,"figures")
if(!dir.exists(figdir)) dir.create(figdir, recursive = T)

################################################################################

# Number of replicates
Iter = 1000

# Assumed design prevalence (district level monthly rate)
Dprev = 1/(100000/12)

buffer_km = 5

# Calculate rolling risk and sensitivity per month during period of interest
tp_m <- c(seq(ymd(cc$start),ymd(cc$end), by = "month"))

# Define matrix of district sensitivity estimates for each time point
sens_cache <- rep(1, length = nrow(shape2))

lapply(tp_m,
       function(tdy){
         run_FFI(tdy,
                 data = data,
                 afp = afp,
                 es = es,
                 R = R,
                 type = type,
                 clust = clust,
                 eg = eg,
                 sens_cache = sens_cache,
                 shapes = shape2,
                 poprast = poprast,
                 proj_local = proj_local,
                 buffer_km = buffer_km,
                 return_catchment_shapes = F,
                 Dprev = Dprev,
                 Iter = Iter) -> out
         # Include to iteratively adjust risk by previous sensitivity
         # sens_cache <<- out$sens$sens_df$SSe_AE
         return(out)
       }) -> out_all

out_all <- setNames(out_all, tp_m)

################################################################################

tp <- sapply(out_all, function(x) as.character(x$sens$sens_period[1]))
tp_f <- format(ymd(tp), "%b-%Y")

# Extract district data for each run
data_all <- lapply(out_all, function(x) x$sens$data) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f))
afp_all <- lapply(out_all, function(x) x$risk$afp) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f))
es_all <- lapply(out_all, function(x) x$risk$es) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f))
sens_all <- lapply(out_all, function(x) x$sens$sens_df) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f))

# Extract ES activity for each run
sens_n_es <- lapply(out_all,
                    function(x) data.frame(Sites = x$sens$n_sites,
                                           Samples = x$sens$n_samples)) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f),
         samples_per_site = Samples/Sites)

# NPEV detection
es_npev_all <- lapply(out_all, function(x) x$risk$es_npev) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f))

# Sensitivity versus coverage
sens_env_wcovg <- lapply(out_all, function(x) as.data.frame(t(x$sens$sens_dist_draws$SeR_ENVi))) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f)) %>%
  bind_cols(select(data_all, guid, es_coverage)) %>%
  select(tdy, tdy_f, guid, es_coverage, everything()) %>%
  pivot_longer(-tdy:-es_coverage, values_to = "SSe_ENV")

# Extract Pneg draws for ENV to recalculate sensitivity in districts with coverage
lapply(out_all,
           function(x) as_tibble(t(x$sens$sens_dist_draws$Pneg_ENVi))) %>%
      setNames(tp) %>%
      bind_rows(.id = "tdy") %>%
      mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                            levels = tp_f)) %>%
      select(tdy, tdy_f, everything()) -> pneg_ENV

# Extract pneg estimates, filter to only LGAs with >X% ES coverage, and recalculate sensitivity overall
perc_covg <- 0.1
data_all %>%
  select(es_coverage) %>%
  bind_cols(pneg_ENV) %>%
  filter(es_coverage > perc_covg) %>%
  pivot_longer(-1:-4, values_to = "Pneg_ENVi", names_to = "Iter") %>%
  group_by(tdy, tdy_f,Iter) %>%
  summarise(SSe_ENV = 1-prod(Pneg_ENVi)) %>%
  ungroup() %>%
  mutate(cat = paste0(">",perc_covg*100,"% coverage")) %>%
  select(-Iter) -> sens_ENV_subset

sens_all %>%
  select(tdy, tdy_f, SSe_ENV) %>%
  mutate(cat = "All") %>%
  bind_rows(sens_ENV_subset) -> sens_env_bycovg

sens_n_afp <- lapply(out_all, function(x) x$sens$data) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  group_by(tdy) %>%
  summarise(across(c(pop_u15_polisdenom,total_pop,afp_n:afp_npev), sum, na.rm = T)) %>%
  ungroup() %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f),
         afp_r = afp_n*1e5/pop_u15_polisdenom)

sens_n_afp_dist <- lapply(out_all, function(x) x$sens$data) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f),
         afp_r = afp_n*1e5/pop_u15_polisdenom)

# Stool adequacy
afp_not_all <- lapply(out_all, function(x) x$risk$afp_not) %>%
  setNames(tp) %>%
  bind_rows(.id = "tdy") %>%
  mutate(tdy_f = factor(format(ymd(tdy),"%b-%Y"),
                        levels = tp_f))

# Save outputs
saveRDS(data_all, here(outdir, "data_all.rds"))
saveRDS(afp_all, here(outdir, "afp_all.rds"))
saveRDS(es_all, here(outdir, "es_all.rds"))
saveRDS(sens_all, here(outdir, "sens_all.rds"))
saveRDS(sens_n_es, here(outdir, "sens_n_es.rds"))
saveRDS(es_npev_all, here(outdir, "es_npev_all.rds"))
saveRDS(sens_env_wcovg, here(outdir, "sens_env_wcovg.rds"))
saveRDS(sens_env_bycovg, here(outdir, "sens_env_bycovg.rds"))
saveRDS(sens_n_afp, here(outdir, "sens_n_afp.rds"))
saveRDS(afp_not_all, here(outdir, "afp_not_all.rds"))
saveRDS(sens_n_afp_dist, here(outdir, "sens_n_afp_dist.rds"))

################################################################################
# Plot maps of each time point
# Prior risk at each time point

lapply(out_all, plot_risk)

# Map samples collected
lapply(out_all, map_samps_env)
# Map estimated sensitivity
lapply(out_all, map_sens_env)

# Map AFP reported
lapply(out_all, map_samps_afp)
# Map estimated sensitivity
lapply(out_all, map_sens_afp)

################################################################################
################################################################################
