################################################################################
################################################################################
# Function to calculate district-specific risk and adjusted risk for given admin2
# data, observed cases/ES samples, and time period
################################################################################

calc_risk <- function(data, 
                      es, 
                      afp,
                      to,
                      from = NULL,
                      retro = NULL,
                      type = "WPV1",
                      clust = NA,
                      eg = NA,
                      prop_u15 = 0.45){
  
# ---------------------------------------------------------------------------- #

# 1) Set up data
  
# A. Subset to time period of interest
  
## Define time period
  to <- ymd(to)
  
  if(!is.null(retro)){
    from = to%m-%months(retro)
  }else if(!is.null(from)){
    from = ymd(from)
  }else{
    return("Provide a value for either <from> or <retro>")
  }

## Subset data
  
  # Subset district data to the year of specified end date 
  # Want one observation per district with an average population/immunity
  data_sub <- data %>% 
    filter(period == date_to_period(floor_date(to, "year"))) %>% 
      select(guid, adm1_name, adm2_name,
             period, population_u5, pop_u15_polisdenom, total_pop, serotype, immunity)
  
  es_sub <- es %>% 
    filter(between(collection_date,from,to)) 

  afp_sub <- afp %>% 
    filter(between(date,from,to)) 
  
  # browser()
## Flag relevant samples/cases 
  if(!is.na(clust)){
    es_sub$flag <- replace_na(grepl(clust, es_sub$virus_cluster_s), FALSE)
    afp_sub$flag <- replace_na(grepl(clust, afp_sub$virus_clusters), FALSE)
  }else if(!is.na(eg)){
    es_sub$flag <- replace_na(grepl(eg, es_sub$emergence_group_s), FALSE)
    afp_sub$flag <- replace_na(grepl(eg, afp_sub$emergence_group_s), FALSE)
  }else{
    es_sub$flag <- replace_na(es_sub$final_class == type, FALSE)
    afp_sub$flag <- replace_na(afp_sub$final_class == type, FALSE)
  }

# B. Aggregate samples/cases

## Summarise AFP cases by district/month: total cases, total positive for WPV1, first/last onset dates for WPV1
  tab_afp <- afp_sub %>% 
    filter(!is.na(guid)) %>% 
    group_by(guid) %>% 
    summarise(afp_n = n(),
              afp_pos=sum(flag, na.rm = T),
              afp_npev=sum(grepl("NPEV",virus_type_s), na.rm = T),
              donset_first_pos=as_date(ifelse(afp_pos!=0,
                                              min(date[flag]),
                                              NA)),
              donset_last_pos=as_date(ifelse(afp_pos!=0,
                                             max(date[flag]),
                                             NA))) %>% 
    ungroup()
  
## Summarise ES by district/month: total samples, total positive for WPV, first/last sample dates
  tab_es <- es_sub %>% 
    filter(!is.na(guid)) %>% 
    group_by(guid) %>%
    summarise(n_sites=n_distinct(site_id),
              es_n=n(),
              es_pos=sum(flag,na.rm=T),
              es_npev=sum(npev == "Yes", na.rm = T),
              collect_first_pos=as_date(ifelse(es_pos!=0,
                                                   min(collection_date[flag]),
                                                   NA)),
              collect_last_pos=as_date(ifelse(es_pos!=0,
                                                  max(collection_date[flag]),
                                                  NA))) %>% 
    ungroup()
  
# C. Add cases/ES to district data
  
  data_sub <- left_join(data_sub,
                        tab_afp)
  
  data_sub <- left_join(data_sub,
                        tab_es)
  
## Set missing counts to 0
  data_sub <- data_sub %>% 
    mutate(across(c("afp_n","afp_pos","afp_npev","n_sites","es_n","es_pos","es_npev"),
                  replace_na, 0))
  
  # browser()
  
# ---------------------------------------------------------------------------- #
# 2) Summarise other sensitivity indicators 
  #  (AFP reporting rate, stool adequacy, ES sampling frequency, Site sampling quality)
  
  # browser()
  
  # Calculate total reporting rate of AFP cases per 100,000 U15 population
  # Compare to target of 2 per 100,000
  # Scale notification parameter by observed rate compared to target: 
  # + Multiplier = 1 for rate >= 2/100,000
  # + Multiplier = observed/target for rate < 2/100,000
  
  # Currently only have district u15 pops for Nigeria. If missing, assume a fixed proportion
  # if (all(is.na(data_sub$pop_u15))){
  #   data_sub <- mutate(data_sub, pop_u15 = prop_u15*total_pop)
  # }
  
  afp_sub %>% 
    group_by(guid) %>% 
    summarise(n_case = n(),
              n_adq = sum(stool_adequacy == "Yes"),
              r_adq = n_adq/n_case,
              r_adq_lwr = Hmisc::binconf(n_adq,n_case)[2],
              r_adq_upr = Hmisc::binconf(n_adq,n_case)[3]) %>% 
    ungroup() %>% 
    right_join(select(data_sub, guid, pop_u15_polisdenom, total_pop)) %>% 
    mutate(n_case = replace_na(n_case, 0),
           r_afp = n_case/pop_u15_polisdenom,
           r_afp_lwr = epitools::pois.exact(n_case,pop_u15_polisdenom)$lower,
           r_afp_upr = epitools::pois.exact(n_case,pop_u15_polisdenom)$upper) -> afp_not
  
  # Sampling and NPEV detection per district
  # Only past *3* months, because NPEV detection is quite seasonal??
  es_sub %>%    
    # filter(between(collection_date, to%m-%months(3), to)) %>% 
    group_by(guid) %>% 
    summarise(n_samp = n(),
              n_site = n_distinct(site_id),
              n_npev = sum(npev == "Yes", na.rm = T),
              r_npev = n_npev/n_samp,
              r_npev_lwr = Hmisc::binconf(n_npev,n_samp)[2],
              r_npev_upr = Hmisc::binconf(n_npev,n_samp)[3]
              ) %>% 
    ungroup() %>% 
    right_join(select(data_sub, guid)) %>% 
    arrange(guid) -> es_npev
  
# ---------------------------------------------------------------------------- #
# 3) Calculate risk

# Note that this is done much better in the regression modelling, and I should 
# really replicate that (?)

  # Default 0 risk
  data_sub$risk_est <- 0

  # Loop through districts
  for(i in 1:dim(data_sub)[1]){
    
    # Cases + ES positives in j, plus radiation
    pow <- (data_sub$afp_pos + data_sub$es_pos)*R[i,]

    # Risk given immunity
    data_sub$risk_est[i] <- 1 - data_sub$immunity[i]^sum(pow)
    
  }
  
  # Log and scale risk value between 0 and 1 for clearer plotting
  # data_sub$risk_est_log <- log(data_sub$risk_est+1E-5)
  # data_sub$risk_est_sc <- replace_na(scale_val(data_sub$risk_est_log), 0)

# ---------------------------------------------------------------------------- #
# 4) Define output
  
  out <- list(data = data_sub, 
              afp = afp_sub,
              es = es_sub,
              risk_period = c(from,to), 
              afp_not = afp_not,
              es_npev = es_npev,
              type = type,
              clust = clust)
  return(out)
  
}


################################################################################
################################################################################