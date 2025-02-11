################################################################################
################################################################################
# Function to estimate sensitivity of AFP and ENV surveillance based on given
# process probabilities, time period of interest and assumed design prevalence
# - include ENV within specific locations
# - include variable importation rate as specified in previous analysis
# - include clustering of infection by having region and unit design prevalence
#
# Questions:
# + Some sort of monthly average coverage rather than just coverage from those 
#   sites sampled this month? 
# + Average over what period?
################################################################################

est_sens <- function(out_risk,
                     es,
                     AFPvals,
                     ESvals,
                     from,
                     to,
                     Iter = 1000, # Iterations of each parameter to incorporate uncertainty
                     Dprev = 1/(100000/12), # Overall design prevalence
                     shapes,
                     poprast,
                     proj_local,
                     buffer_km = 5,
                     return_catchment_shapes = F,
                     equal_risk = F,
                     retro = T){

  from <- ymd(from)
  to <- ymd(to)

  data <- out_risk$data
  numDist <- dim(data)[1]
  
  afp_not <- out_risk$afp_not
  es_npev <- out_risk$es_npev
  
  # browser()
  
  # U15 population proportion in each node/district
  nR <- data$pop_u15_polisdenom
  PrP <- nR/sum(nR)         
  
  # Design prevalence
  Dprev_region <- rep(1/numDist, Iter)   # should detect 1 region in the total population (regional design prevalence)                             
  
  Dprev_unit <- rep(Dprev,Iter) 

  # ---------------------------------------------------------------------------- #
  # Function to draw from beta distribution defined by quantiles  

  get_beta_draws <- function(x, p = 0.95, Iter = 1000){
    
    if(!any(is.na(x))){
      tmp <- betaExpert(best=x[1],
                        lower=x[2],
                        upper=x[3],p=p)
    }
    draws <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta) 
    return(draws)
  }
  
  # Helper function as betaExpert can't handle exact 0/1s
  avoid_01 <- function(x){
    x[x <= 0] <- 0.00001
    x[round(x,5) >= 1] <- 0.99999
    return(x)
  }
  
  # Adjust observed rates relative to a specified target
  target_scale <- function(x, tgt){
    x <- ifelse(x < tgt, x/tgt, 1)
  }
  
  ################################################################################
  # STEP 1: go through each of detection sensitivities
  
  # ---------------------------------------------------------------------------- #
  # (A). AFP surveillance
  
  # Incorporating uncertainty in assumptions of AFP surveillance sensitivity
  
  # Infer beta distribution parameters from point estimate + interval for each detection step:
  # - Rate of AFP
  # - P[notified]
  # - P[adequate stool sample]
  # - P[test positive]
  # - sensitivity of the test
  
  # Draw from distributions:

  #----------------------------------------------------------------------#
  # P_AFPclinical_wild - Probability of developing paralysis if infected #
  #----------------------------------------------------------------------#
  
  tmp <- betaExpert(best=1/AFPvals$EstMean[which(AFPvals$SurveillanceNode=="AFPcase")],
                    lower=1/AFPvals$EstLwr[which(AFPvals$SurveillanceNode=="AFPcase")],
                    upper=1/AFPvals$EstUpr[which(AFPvals$SurveillanceNode=="AFPcase")],p=0.95)
  P_AFPclinical_wild <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # CLINICAL case (<1%) ...  # hist(Naf)
  # signif(quantile(P_AFPclinical_wild,probs=c(0.025,0.5,0.975)),3)[c(2,1,3)] 
  
  #------------------------------------------------------#
  # P_AFPnot - Probability of an AFP case being notified #
  #------------------------------------------------------#
  
  # General estimate
  tmp <- betaExpert(best=AFPvals$EstMean[which(AFPvals$SurveillanceNode=="AFPnotified")],
                    lower=AFPvals$EstLwr[which(AFPvals$SurveillanceNode=="AFPnotified")],
                    upper=AFPvals$EstUpr[which(AFPvals$SurveillanceNode=="AFPnotified")],p=0.95)
  P_AFPnot <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # Notification *rate* (this is still "uncertain")
  
  # Scale by previous observed AFP rate per district
  # Minority of districts reporting < 2/100,000 per year 
  
  # Fill in districts without obs with an average across other districts
  afp_not$r_afp[is.na(afp_not$r_afp)] <- median(afp_not$r_afp, na.rm = T)
  afp_not$r_afp_lwr[is.na(afp_not$r_afp_lwr)] <- median(afp_not$r_afp_lwr, na.rm = T)
  afp_not$r_afp_upr[is.na(afp_not$r_afp_upr)] <- median(afp_not$r_afp_upr, na.rm = T)
  
  # Define scaled value according to specified target rate
  target_afp = 2/1e5
  afp_not %>% 
    mutate(across(r_afp:r_afp_upr, 
                  function(x) avoid_01(target_scale(x, tgt = target_afp)))) -> afp_not_scaled
  
  # Add beta variation for each district
  P_AFPnoti <- apply(select(afp_not_scaled,r_afp,r_afp_lwr,r_afp_upr),1,get_beta_draws)
  
  P_AFPnoti <- P_AFPnoti*P_AFPnot
  
  #-------------------------------------------------------#
  # P_AFPstooli - Probability of adequate stool collected #
  #-------------------------------------------------------#

  # Observed during risk period
  afp_not %>% 
    pivot_longer(r_adq:r_adq_upr) %>% 
    ggplot(aes(value, fill = name)) +
    facet_wrap(~name) +
    guides(fill = "none") +
    geom_histogram(alpha = 0.5, position = "identity") +
    labs(x = "Proportion", y = "Freq",
         title = "Proportion of adequate stool samples by district",
         subtitle = "Mean and 95% binomial confidence limits",
         caption = paste(out_risk$risk_period, collapse = " - ")) -> p_stool_vals
  
  # Fill in districts without obs with average across other districts
  afp_not$r_adq[is.na(afp_not$r_adq)] <- median(afp_not$r_adq, na.rm = T)
  afp_not$r_adq_lwr[is.na(afp_not$r_adq_lwr)] <- median(afp_not$r_adq_lwr, na.rm = T)
  afp_not$r_adq_upr[is.na(afp_not$r_adq_upr)] <- median(afp_not$r_adq_upr, na.rm = T)
  
  # Add beta variation for each district  
  afp_not %>% 
    mutate(across(r_adq:r_adq_upr, 
                  function(x) avoid_01(x))) -> afp_stool_scaled
  P_AFPstooli <- apply(select(afp_stool_scaled,r_adq,r_adq_lwr,r_adq_upr),1,get_beta_draws)
  
  #------------------------------------------------------------#
  # P_AFPtest - Probability of positive result given infection #
  #------------------------------------------------------------#
  
  tmp <- betaExpert(best=AFPvals$EstMean[which(AFPvals$SurveillanceNode=="AFPtest")],
                    lower=AFPvals$EstLwr[which(AFPvals$SurveillanceNode=="AFPtest")],
                    upper=AFPvals$EstUpr[which(AFPvals$SurveillanceNode=="AFPtest")],p=0.95)
  P_AFPtest <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   
  
  #-------------------------------------------------------#
  
  # Plot sampled parameters:
  P_AFPnot <- P_AFPnoti %>% 
    as_tibble() %>% 
    pivot_longer(everything()) %>% 
    mutate(name = "P_AFPnoti")
  P_AFPstool <- P_AFPstooli %>% 
    as_tibble() %>% 
    pivot_longer(everything()) %>% 
    mutate(name = "P_AFPstooli") 

  bind_cols(P_AFPclinical_wild,
            P_AFPtest) %>% 
    setNames(c("P_AFPclinical_wild",
               "P_AFPtest")) %>% 
    pivot_longer(everything()) %>%
    bind_rows(P_AFPstool) %>% 
    bind_rows(P_AFPnot) %>% 
    mutate(name = factor(name, levels = c("P_AFPclinical_wild","P_AFPnoti","P_AFPstooli","P_AFPtest"))) %>% 
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(~name, scales = "free_y") -> plots_afp
  
  # ---------------------------------------------------------------------------- #
  # (B). ENV surveillance
  
  # Sampling activity within this time period
  es_sub <- es %>% 
    filter(between(collection_date, from, to))

  #--------------------------------------------------------#
  # P_ENVi - Probability of being within ES catchment area #
  #--------------------------------------------------------#
  
  covg_out <- get_es_coverage(data = data,
                              es_linelist = es_sub,
                              shapes = dplyr::select(shapes, guid, total_pop,
                                                     center_lon, center_lat),
                              poprast = poprast,
                              proj = proj_local,
                              buffer_km = buffer_km,
                              return_catchment_shapes = return_catchment_shapes)
  P_ENVi <- covg_out$data$es_coverage
  
  # Also add ES catchment variables to data
  data <- covg_out$data %>%
      right_join(dplyr::select(shapes, guid, total_pop) %>% st_drop_geometry())
  
  # NOTE: Sensitivity analysis varying catchment radius demonstrates impact of radius choice
    
  # browser()
  # Add beta variation for each district
  tmp <- data.frame(est = P_ENVi, lwr = P_ENVi*0.85, upr = P_ENVi*1.15) %>% 
    apply(2,avoid_01)
  P_ENVi <- apply(tmp,1,get_beta_draws)
  
  #------------------------------------------------------#
  # P_ENVsampi - Probability of sampling during shedding #
  #------------------------------------------------------#
  
  # If the analysis is retrospective, we can estimate based on each month's observed
  # sampling. 
  # If it is prospective (i.e. we want to estimate sensitivity as of now and assume that persists for future FFI),
  # might be better to use an average frequency across the past 12?
  
  if (retro == F){
    # Average monthly sampling frequency across all districts during prior 12m
    freq_est <- mean(es_freq$avg_freq_mth, na.rm = T)
    
    # -> Probability of sampling during shedding OVERALL
    P_ENVsamp <- get_p_envsamp(freq_est)
    
    # Average monthly sampling frequency per district during prior 12m
    freq_esti <- es_freq$avg_freq_mth
    
    table(is.na(freq_esti))
    
    # Note: this doesn't work if a site is activated this month, but has no data from
    # prior months with which to estimate frequency
    # -> if new site assume once a month? conservative
    freq_esti[is.na(freq_esti)] <- 1 
    
  }else{
    # Number of samples collected per district THIS month
    es_sub %>% 
      group_by(guid) %>% 
      summarise(n = n()) %>% 
      ungroup() %>% 
      right_join(select(data, guid)) %>% 
      mutate(n = replace_na(n, 0)) %>%
      arrange(guid) %>% pull(n) -> freq_esti
  }

  # -> Probability of sampling during shedding per district
  P_ENVsampi <- get_p_envsamp(freq_esti)
  
  #--------------------------------------------------------#
  # P_ENVqual - Probability of sufficient quality sampling #
  #--------------------------------------------------------#

  # Fill in districts without obs with an average across other districts
  es_npev$r_npev[is.na(es_npev$r_npev)] <- median(es_npev$r_npev, na.rm = T)
  es_npev$r_npev_lwr[is.na(es_npev$r_npev_lwr)] <- median(es_npev$r_npev_lwr, na.rm = T)
  es_npev$r_npev_upr[is.na(es_npev$r_npev_upr)] <- median(es_npev$r_npev_upr, na.rm = T)

  # Define scaled value according to specified target rate
  target_npev = 0.5
  es_npev %>% 
    mutate(across(r_npev:r_npev_upr, 
                  function(x) avoid_01(target_scale(x, tgt = target_npev)))) -> es_npev_scaled
  
  # Add beta variation for each district
  P_ENVquali <- apply(select(es_npev_scaled,r_npev,r_npev_lwr,r_npev_upr),1,get_beta_draws)

  #------------------------------------------------------------#
  # P_ENVtest - Probability of positive result given infection #
  #------------------------------------------------------------#
  
  tmp <- betaExpert(best=ESvals$EstMean[which(ESvals$SurveillanceNode=="EStest")],
                    lower=ESvals$EstLwr[which(ESvals$SurveillanceNode=="EStest")],
                    upper=ESvals$EstUpr[which(ESvals$SurveillanceNode=="EStest")],p=0.95)
  P_ENVtest <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # ES test sens (this might be low because of CSF)
  
  #--------------------------------------------------------#
  
  # Plot sampled parameters:
  data.frame(name = "P_ENVi",value = as.vector(P_ENVi)) %>% 
    bind_rows(data.frame(name = "P_ENVsampi",value = as.vector(P_ENVsampi))) %>% 
    bind_rows(data.frame(name = "P_ENVquali",value = as.vector(P_ENVquali))) %>% 
    bind_rows(data.frame(name = "P_ENVtest", value = P_ENVtest)) %>% 
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(~name, scales = "free_y") +
    labs(caption = "Including zero coverage districts")-> plots_env1
  
  data.frame(name = "P_ENVi",value = as.vector(P_ENVi)) %>% 
    bind_rows(data.frame(name = "P_ENVsampi",value = as.vector(P_ENVsampi))) %>% 
    bind_rows(data.frame(name = "P_ENVquali",value = as.vector(P_ENVquali))) %>% 
    bind_rows(data.frame(name = "P_ENVtest", value = P_ENVtest)) %>%
    filter(value != 0) %>% 
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(~name, scales = "free_y") +
    labs(caption = "Excluding zero coverage districts") -> plots_env2
  
  plots_env <- plots_env1 / plots_env2
  plots_env
  
  ################################################################################
  # STEP 2: based on above, calculate Unit Sensitivity (USe) and Surveillance 
  # Sensitivity (SSe) for 1 month.
  
  # Accounting for DIFFERENTIAL RISK
  
  if(equal_risk == T){
    data$risk_est <- 1
    data$risk_est_log <- 0
    data$risk_est_sc <- 1
  }else{
    # Log and scale estimated risk
    # Specify scaled risk as 0 if original risk was 0
    data$risk_est_log <- na_if(log(data$risk_est),-Inf)
    data$risk_est_sc <- replace_na(scale_val(data$risk_est_log), 0)
  }

  tmp <- data$risk_est_sc+1e-5      #+1e-10 #(plus a bit so min risk >0.00) #  a log scale and not min=1

  # browser()
  
  # Define relative risk between districts, setting lowest risk as the reference level
  RR <- tmp/min(tmp)
  hist(RR) 
  
  # Solve for adjusted risk (AR) based on population proportions
  data_AR <- data.frame(RR=RR, PrP=PrP)
  AR <-  ARfunct2(data_AR)
  
  # Repeat for population covered by ES
  # Assume same proportion of total U15s and total pop covered by ENV 
  # nR_ENV <- data$pop_u15*P_ENVi
  # PrP_ENV <- nR_ENV/sum(nR_ENV)
  
  # data_AR_ENV <- data.frame(RR=RR, PrP=PrP_ENV)
  # AR_ENV <-  ARfunct2(data_AR_ENV)
  
  # browser()
  
  # Check constraints
  if(sum(PrP*AR) != 1.00){ #| sum(PrP_ENV*AR_ENV) != 1
    print("Warning: Adjusted risks not equal to 1")
  }

  # Add AR to dataset
  data$AR <- AR
  
  # ---------------------------------------------------------------------------- #
  # Sensitivity across scenario trees
  
  # "Unit" sensitivity of AFP: Probability of detection given infection
  # - Probability of clinical symptoms (fixed)
  # - Probability of being notified (fixed)
  # - Probability of having adequate stool collected (average proportion + CI for risk period, over districts)
  # - Probability of test positive given infection (fixed)
  
  sens_AFPi <- P_AFPnoti*P_AFPstooli*matrix(rep(P_AFPclinical_wild*P_AFPtest, ncol(P_AFPstooli)),
                                            nrow = Iter)
  # Summarise across draws
  signif(quantile(sens_AFPi,
                  probs = c(0.5,0.025,0.975)),3) 
  
  hist(sens_AFPi,
       probability=T,
       xlab = "Unit sensitivity of Detecting 1 Infection",
       main="AFP Surveillance")
  
  # "Unit" sensitivity of ENV:
  # - Probability of being in catchment area (per district)
  # - Probability of sample collection (fixed)
  # - Probability test positive given infection (fixed)
  
  ### Previously no variation by district, only from iterations -> 1000 vector
  ### Now also coverage depends on district so sampling variation applied to each district separately -> nDist*1000 matrix
  ### Repeat each district covg for each iteration, and each test sensitivity iteration for each district
  # sens_ENVi <- matrix(rep(P_ENVi, each = Iter),nrow = Iter)*P_ENVsampi*P_ENVquali*matrix(rep(P_ENVtest,length(P_ENVi)),nrow = Iter)
  sens_ENVi <- P_ENVi*P_ENVsampi*P_ENVquali*matrix(rep(P_ENVtest,ncol(P_ENVi)),nrow = Iter)
  
  # Overall sensitivity (incl variation bw districts)
  signif(quantile(sens_ENVi,
                  probs = c(0.5,0.025,0.975)),3)
  hist(sens_ENVi,
       probability=T,
       xlab = "Unit sensitivity of Detecting 1 Infection",
       main="ES Surveillance")
  
  # Sensitivity within districts with coverage
  signif(quantile(sens_ENVi[,which(colSums(P_ENVi) > 0)],
                  probs = c(0.5,0.025,0.975)),3)
  hist(sens_ENVi[,which(colSums(P_ENVi) > 0)],
       probability=T,
       xlab = "Unit sensitivity of Detecting 1 Infection, given ES coverage",
       main="ES Surveillance")
  
  # Combined AFP and ENV (wild) - Non-overlapping so just add together
  # signif(quantile(sens_AFPi + sens_ENVi, probs = c(0.5,0.025,0.975)),3)
  # 
  # hist(sens_AFPi + sens_ENVi,
  #      probability=T,
  #      xlab = "Unit sensitivity of Detecting 1 Infection",
  #      main="AFP + ES Surveillance")
  
  # ---------------------------------------------------------------------------- #
  # Calculate region-specific sensitivity
  
  # Combine surveillance with the effective probability of infection (EPI):
  # EPI = AR*Ph (if only thinking about Ph - 1/100,000 in UK % 12? for per month?)
  # EPI = AR*Ph*Pr  (detecting 1 infected region with at least 1/100000 per region)
  # SeR_AFPi = Sensitivity of surveillance (AFP) per district for 1 month of data with negatives
  # Pneg_AFPi = Probability that a district will have a negative result
  #             ... 1-Pr(1+ pos result)    accounting for risk (This isn't that useful for district lvel analysis)
  
  
  #######
  # AFP #
  #######
  
  # Setup outputs: matrix of iterations (rows) and districts (cols)
  SeR_AFPi <- Pneg_AFPi <- matrix(0,nrow=Iter, ncol=numDist)
  
  for(i in 1:numDist){
    
    # District sensitivity from unit design prevalence, overall AFP sensitivity, and population size:
    # (varying risk of infection by age group would be added here, with a product of the inner (1 - ...) over risk groups)
    SeR_AFPi[,i] <- 1 - (1 - (Dprev_unit*P_AFPclinical_wild*P_AFPnoti[,i]*P_AFPstooli[,i]*P_AFPtest))^nR[i]  
    
    # Probability per district of negative result from region design prevalence, region sensitivity and adjusted risk:
    Pneg_AFPi[,i] <- 1 - (Dprev_region*SeR_AFPi[,i]*AR[i])
    
  }

  # Average sensitivity per district (accounting for EPI)
  tmp <- apply(SeR_AFPi,2,mean) 
  
  # Add to dataset
  data$SeR_AFPi <- tmp

  # Sensitivity across the entire system (1 - product of region negative result probabilities)
  # - Per draw
  SSCSe_AFP <- rep(NA,Iter)
  for(j in 1:Iter){
    SSCSe_AFP[j] <- 1 - prod(Pneg_AFPi[j,])
  }
  
  #######
  # ENV #
  #######
  # - assuming different ES coverage per district dependent on time period
  
  # Setup outputs
  SeR_ENVi <- Pneg_ENVi <- matrix(0,nrow=Iter,ncol=numDist)  
  
  # Power is the number of units for which we have negative outcomes, i.e. the number of individuals in ES catchment
  for(i in 1:numDist){
    # accounting for different coverage (P_ENV) in each district
    SeR_ENVi[,i] <- 1 - (1 - (Dprev_unit*P_ENVi[,i]*P_ENVsampi[,i]*P_ENVtest))^nR[i]
    Pneg_ENVi[,i] <- 1 - (Dprev_region*SeR_ENVi[,i]*AR[i])
  }

  # Average sensitivity per district
  tmp <- apply(SeR_ENVi,2,mean)
  # NOTE: sens = 0/Pneg = 1 in districts with no ES site. 
  
  # Add to dataset
  data$SeR_ESi <- tmp
  
  # ggplot(datdist, aes(es_coverage, SeR_ESi)) +
  #   geom_point(alpha = 0.2) +
  #   geom_smooth() + 
  #   scale_x_continuous(trans = "sqrt") + 
  #   scale_y_continuous(trans = "sqrt")
  # 
  # ggplot(datdist, aes(AR_ENV, SeR_ESi)) +
  #   geom_point(alpha = 0.2) +
  #   geom_smooth() +
  #   scale_y_continuous(trans = "sqrt") 
  
  # Sensitivity across the entire system (1 - product of region negative result probabilities)
  # - Per draw
  SSCSe_ENV <- rep(NA,Iter)
  for(j in 1:Iter){
    SSCSe_ENV[j] <- 1 - prod(Pneg_ENVi[j,])
  }

  ############
  # Combined #
  ############
  
  # Combining the surveillance together....
  SSe_AFP <- 1 - ((1-SSCSe_AFP))
  
  SSe_ENV <- 1 - ((1-SSCSe_ENV))
  
  # Assumed independent so multiply
  SSe_AE <- 1 - ((1-SSCSe_AFP)*(1-SSCSe_ENV))
  
  ##############################################################################
  # Define output
  
  sens_df <- data.frame(SSe_AFP = SSe_AFP, 
                        SSe_ENV = SSe_ENV, 
                        SSe_AE = SSe_AE)
  
  sens_dist_draws <- list(SeR_AFPi = SeR_AFPi,
                          Pneg_AFPi = Pneg_AFPi,
                          SeR_ENVi = SeR_ENVi, 
                          Pneg_ENVi = Pneg_ENVi)

  out <- list(sens_df = sens_df, 
              sens_dist_draws = sens_dist_draws,
              data = data,
              sens_period = c(from,to),
              n_sites = covg_out$n_sites,
              n_samples = covg_out$n_samples,
              # afp_not_draws = P_AFPstooli,
              # es_sampfreq_draws = P_ENVsampi,
              plots_afp = plots_afp,
              plots_env = plots_env,
              catchment_shapes = covg_out$catchment_shapes)
  return(out)

}


# est_sens_all <- lapply(tp, est_sens)
# saveRDS(est_sens_all, here(outdir,"est_sens_all.rds"))

################################################################################
################################################################################