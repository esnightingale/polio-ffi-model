################################################################################
################################################################################
# Function for entire workflow
################################################################################


run_FFI <- function(tdy,
                    retro = 12,
                    type = "WPV1",
                    clust = NULL,
                    eg = NULL,
                    data,
                    afp,
                    es,
                    R,
                    shapes,
                    poprast,
                    proj_local = proj_local,
                    buffer_km = 5,
                    return_catchment_shapes = F,
                    equal_risk = F,
                    Dprev = 1/(100000/12),
                    Iter = 1000,
                    sens_cache = sens_cache,
                    est_FFI = F,
                    M,
                    elim_criterion = 0.95,
                    seed = 1111){

  # -----------------------------------------------------------------------------#
  # SOURCE SUB-FUNCTIONS

  source("R/utils/calculate_risk_fcn.R")
  source("R/utils/estimate_sens_fcn.R")
  source("R/utils/estimate_elimination_fcn.R")

  source("R/utils/get_es_coverage.R")
  source("R/utils/mapping_fcns.R")
  source("R/utils/helper_fcns.R")

  # Today
  tdy = ymd(tdy)

  # Set seed
  set.seed(seed)

  # Check M
  if(est_FFI & is.null(M)){
    return("Please provide a value for M")
  }

  # -----------------------------------------------------------------------------#
  # CALCULATE RISK

  # Define period across which to calculate risk
  # - Either define "retro" months prior to this month
  # - Or define particular "from" date

  risk_to = floor_date(tdy-1, "month")
  risk_from = risk_to%m-%months(retro)

  # Alternatively:
  # if(is.null(retro)){
  #   risk_from = ymd(risk_from)
  # }

  print(paste("Risk calculated from",risk_from, "to",risk_to))

  # if (!first_flag){}

  # Calculate risk per district based on observations within the specified period
  calc_risk(data,
            es,
            afp,
            from = risk_from,
            to = risk_to,
            retro = retro,
            type = type,
            clust = clust,
            eg = eg,
            sens_last = sens_cache) -> out_risk

  # browser()

  # -----------------------------------------------------------------------------#
  # ESTIMATE SENSITIVITY

  # Define period across which to estimate sensitivity
  sens_from = floor_date(tdy-1,"month")
  sens_to = tdy

  print(paste("Sensitivity calculated from",sens_from,"to",sens_to))

  # Estimate sensitivity based on sampling activity during *current* period
  est_sens(out_risk = out_risk,  # District population and risk data (from previous step)
           es = es,               # *Un-subsetted* ES samples
           AFPvals = AFPvals,
           ESvals = ESvals,
           from = sens_from,
           to = sens_to,
           Iter = Iter,
           Dprev = Dprev, # Overall design prevalence per month
           shapes = shape2,
           poprast = poprast,
           proj_local = proj_local,
           buffer_km = buffer_km,
           return_catchment_shapes = return_catchment_shapes,
           equal_risk = equal_risk) -> out_sens

  # browser()

  # Plot catchment areas for an example region
  # xlim = c(66.8,67.3)
  # ylim = c(24.7,25.3)
  # lims = st_as_sf(x = tibble(x = xlim, y = ylim),
  #                 coords = c("x","y"),
  #                 crs = st_crs(4326)) %>%
  #   st_transform(crs = proj_local)
  # ggplot() +
  #     geom_sf(data = shape2) +
  #     geom_sf(data = out_sens$catchment_shapes, aes(fill = guid.1)) +
  #     coord_sf(xlim = st_coordinates(lims)[,1], ylim = st_coordinates(lims)[,2]) +
  #     guides(fill = "none")

  # -----------------------------------------------------------------------------#
  # ESTIMATE FREEDOM FROM INFECTION

  # Estimate FFI based on calculated sensitivities for **current** month
  # For M **prospective** years of zero positives
  #
  # Note: It's probably unnecessary repetition to calculate an M-year-ahead
  # FFI trajectory for each month's estimated sensitivity. Usually want to set
  # est_FFI = F and instead calculate a single trajectory over a retrospective
  # period with variable sensitivity via 05_estimate_elimination_trajectory.R.
  # This can then be projected into the future based on last 3m average
  # sensitivities by setting future = T and providing M.

  if(est_FFI == T){

    # Define period across which to estimate FFI
    ffi_from = sens_from #tdy
    ffi_to = ffi_from+years(M)

    print(paste("FFI calculated assuming zero positives from",ffi_from, "to",ffi_to))

    elim_est(out_sens = out_sens,
             tdy = ffi_from,
             future = T,
             M = M,
             elim_criterion = elim_criterion) -> out_elim
    out <- list(risk = out_risk,
                sens = out_sens,
                FFI = out_elim)

  }else{
    out <- list(risk = out_risk,
                sens = out_sens)
  }

  return(out)

}
