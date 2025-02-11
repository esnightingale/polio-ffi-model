elim_est_tv <- function(sens_all,
                     seed = 1111,
                     tv = T,
                     future = F,
                     M = NA,
                     prior_parms = c(30,30), 
                     intro = T,
                     intro_parms = c(lower = 1/5000, 
                                     best = 1/1000, 
                                     upper = 1/500, 
                                     p = 0.95), # default low probability of introduction
                     elim_criterion = 0.95){
  
  set.seed(seed)
  
  # Extracted estimated sensitivities for each run
  sens_all %>% 
    group_by(tdy,tdy_f) %>% 
    mutate(iter = row_number()) %>% 
    ungroup() -> sens_all
  
  sens_all %>% 
    select(tdy,iter,SSe_AFP) %>% 
    pivot_wider(values_from = SSe_AFP, names_from = iter) %>% 
    column_to_rownames("tdy") %>% 
    as.matrix() -> SSe_AFP
  sens_all %>% 
    select(tdy,iter,SSe_ENV) %>% 
    pivot_wider(values_from = SSe_ENV, names_from = iter) %>% 
    column_to_rownames("tdy") %>% 
    as.matrix() -> SSe_ENV
  sens_all %>% 
    select(tdy,iter,SSe_AE) %>%  
    pivot_wider(values_from = SSe_AE, names_from = iter) %>% 
    column_to_rownames("tdy") %>% 
    as.matrix() -> SSe_AE
  
  Iter = ncol(SSe_AE)
  max_time = nrow(SSe_AE)    
  last_time <- ymd(last(rownames(SSe_AFP)))
  
  if(!tv){
    for (i in 2:nrow(SSe_AFP)) {
      SSe_AFP[i, ] <- SSe_AFP[1, ]
    }   
    for (i in 2:nrow(SSe_ENV)) {
      SSe_ENV[i, ] <- SSe_ENV[1, ]
    }
    for (i in 2:nrow(SSe_AE)) {
      SSe_AE[i, ] <- SSe_AE[1, ]
    }
  }
  
  if(future){
    
    if(is.null(M)){print("Provide a value for M")}
    add_time <- M*12
    max_time <- max_time + add_time

    tp <- seq(last_time + months(1),
                last_time + months(add_time),
                by = "month")
    
    # Pull average of last 3m estimated sensitivities
    SSe_AFP_1 <- colMeans(SSe_AFP[(nrow(SSe_AFP)-2):nrow(SSe_AFP),]) %>% 
      matrix(byrow = T, ncol = Iter, nrow = add_time)
    rownames(SSe_AFP_1) <- as.character(tp)
    SSe_ENV_1 <- colMeans(SSe_ENV[(nrow(SSe_ENV)-2):nrow(SSe_ENV),]) %>% 
      matrix(byrow = T, ncol = Iter, nrow = add_time)
    rownames(SSe_ENV_1) <- as.character(tp)
    SSe_AE_1 <- colMeans(SSe_AE[(nrow(SSe_AE)-2):nrow(SSe_AE),]) %>% 
      matrix(byrow = T, ncol = Iter, nrow = add_time)
    rownames(SSe_AE_1) <- as.character(tp)
    
    # Define new matrices including future time points
    SSe_AFP = rbind(SSe_AFP,SSe_AFP_1)
    SSe_ENV = rbind(SSe_ENV,SSe_ENV_1)
    SSe_AE = rbind(SSe_AE, SSe_AE_1)

  }

  # head(SSe_AFP)
  # head(SSe_AE)   
  
  # TEST to see impact of doubling sensitivity
  # }else if(tv){
  #   SSe_AFP = 2*SSe_AFP
  #   SSe_AE = 2*SSe_AE
  # }
  
  # ----------------------------------------------------------------------------#
  # Setup outputs - matrices of time points by iteration
  
  time <- 0:max_time
  
  # SCENARIO A: NO INTRODUCTION RISK
  ## Both 
  Post_freeA <- matrix(0,nrow=length(time),ncol=Iter)   
  Prior_infeA <- matrix(0,nrow=length(time),ncol=Iter)
  
  ## AFP surveillance only
  Post_freeA_AFP <- matrix(0,nrow=length(time),ncol=Iter)   
  Prior_infeA_AFP <- matrix(0,nrow=length(time),ncol=Iter)
  
  ## ENV surveillance only
  Post_freeA_ENV <- matrix(0,nrow=length(time),ncol=Iter)  
  Prior_infeA_ENV <- matrix(0,nrow=length(time),ncol=Iter)
  
  # SCENARIO B: WITH INTRODUCTION RISK
  ## Both
  Post_freeB <- matrix(0,nrow=length(time),ncol=Iter) 
  Prior_infeB <- matrix(0,nrow=length(time),ncol=Iter)
  
  ## AFP Surveillance only
  Post_freeB_AFP <- matrix(0,nrow=length(time),ncol=Iter) 
  Prior_infeB_AFP <- matrix(0,nrow=length(time),ncol=Iter)
  
  ## ENV surveillance only
  Post_freeB_ENV <- matrix(0,nrow=length(time),ncol=Iter)  
  Prior_infeB_ENV <- matrix(0,nrow=length(time),ncol=Iter)
  
  # -------------------------------------------------------------------------- #
  # Probability of introduction
  # - Currently fixed value for all districts
  # - Should be updated to be region-specific
  # - Proportional to population size?
  
  # default : best=1/(1000),lower=1/5000,upper=1/500
  # should change it really...but it's small mostly
  
  PP_intro <-   betaExpert(best = intro_parms["best"],
                           lower = intro_parms["lower"],
                           upper = intro_parms["upper"], p = intro_parms["p"])  # low
  # PP_intro <-   betaExpert(best=1/(100),lower=1/1000,upper=1/10,p=0.95)  # high risk
  #PP_intro <-   betaExpert(best=1/(1000),lower=1/5000,upper=1/500,p=0.95)  # default
  
  # Draws from this distribution
  Pr_intro <- rbeta(n=Iter,shape1=PP_intro$alpha,shape2=PP_intro$beta)
  hist(Pr_intro) 
  
  # -------------------------------------------------------------------------- #
  # Updating FFI per month without observed disease
  
  # initially, we are unsure about freedom - prior is 0.5
  # let's say last case was 6 mths ago
  # if as a country we assume a 50% chance, we need to split it into the 
  # different LAs - otherwise risk is high!
  
  # Histogram of assumed prior distribution for FFI
  hist(rbeta(Iter,prior_parms[1],prior_parms[2]),xlim=c(0,1))
  
  survi <- "wild"  #"vdpv"  
  if(survi=="wild"){
    Prior_infeA[1,] <- Prior_infeB[1,] <- Post_freeA[1,] <- Post_freeB[1,] <- rbeta(Iter,prior_parms[1],prior_parms[2])
    Prior_infeA_AFP[1,] <- Prior_infeB_AFP[1,] <- Post_freeA_AFP[1,] <- Post_freeB_AFP[1,] <- rbeta(Iter,prior_parms[1],prior_parms[2])
    Prior_infeA_ENV[1,] <- Prior_infeB_ENV[1,] <- Post_freeA_ENV[1,] <- Post_freeB_ENV[1,] <- rbeta(Iter,prior_parms[1],prior_parms[2])
  }
  
  # start <- date_to_period(tdy)
  start <- date_to_period(min(sens_all$tdy)) 
  
  for(t in 1:(max_time)){
    
    # timei <- floor(t/12)+start  # time in years
    print(t)
    
    # SCENARIO A: Ignoring introduction risk 
    
    # BOTH surveillance components
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeA[t+1,] <- (1-(1-Post_freeA[t,]))/(1-(1-Post_freeA[t,])*SSe_AE[t,sample(1:Iter,Iter)])
    
    # AFP SURVEILLANCE
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeA_AFP[t+1,] <- (1-(1-Post_freeA_AFP[t,]))/(1-(1-Post_freeA_AFP[t,])*SSe_AFP[t, sample(1:Iter,Iter)]) #
    
    # ENV SURVEILLANCE
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeA_ENV[t+1,] <- (1-(1-Post_freeA_ENV[t,]))/(1-(1-Post_freeA_ENV[t,])*SSe_ENV[t, sample(1:Iter,Iter)]) #
    
    # SCENARIO B: Including introduction risk 
    # Rearranging definition of prior infect
    
    # # BOTH surveillance components
    # Post_freeB[t+1,] <- (1-(Prior_infeB[t,]))/(1-(Prior_infeB[t,])*SSe_AE[t, sample(1:Iter,Iter)])
    # Prior_infeB[t+1,] <- (1-Post_freeB[t+1,]) + Pr_intro[sample(1:Iter,Iter)] - ((1-Post_freeB[t+1,])*Pr_intro[sample(1:Iter,Iter)])
    # 
    # # AFP SURVEILLANCE
    # Post_freeB_AFP[t+1,] <- (1-(Prior_infeB_AFP[t,]))/(1-(Prior_infeB_AFP[t,])*SSe_AFP[t, sample(1:Iter,Iter)])
    # Prior_infeB_AFP[t+1,] <- (1-Post_freeB_AFP[t+1,]) + Pr_intro[sample(1:Iter,Iter)] - ((1-Post_freeB_AFP[t+1,])*Pr_intro[sample(1:Iter,Iter)])
    # 
    # # ENV SURVEILLANCE
    # Post_freeB_ENV[t+1,] <- (1-(Prior_infeB_ENV[t,]))/(1-(Prior_infeB_ENV[t,])*SSe_ENV[t, sample(1:Iter,Iter)])
    # Prior_infeB_ENV[t+1,] <- (1-Post_freeB_ENV[t+1,]) + Pr_intro[sample(1:Iter,Iter)] - ((1-Post_freeB_ENV[t+1,])*Pr_intro[sample(1:Iter,Iter)])
    
    # BOTH surveillance components
    Prior_infeB[t,] <- (1-Post_freeB[t,]) + Pr_intro[sample(1:Iter,Iter)] - ((1-Post_freeB[t,])*Pr_intro[sample(1:Iter,Iter)])
    Post_freeB[t+1,] <- (1-(Prior_infeB[t,]))/(1-(Prior_infeB[t,])*SSe_AE[t, sample(1:Iter,Iter)])

    # AFP SURVEILLANCE
    Prior_infeB_AFP[t,] <- (1-Post_freeB_AFP[t,]) + Pr_intro[sample(1:Iter,Iter)] - ((1-Post_freeB_AFP[t,])*Pr_intro[sample(1:Iter,Iter)])
    Post_freeB_AFP[t+1,] <- (1-(Prior_infeB_AFP[t,]))/(1-(Prior_infeB_AFP[t,])*SSe_AFP[t, sample(1:Iter,Iter)])

    # ENV SURVEILLANCE
    Prior_infeB_ENV[t,] <- (1-Post_freeB_ENV[t,]) + Pr_intro[sample(1:Iter,Iter)] - ((1-Post_freeB_ENV[t,])*Pr_intro[sample(1:Iter,Iter)])
    Post_freeB_ENV[t+1,] <- (1-(Prior_infeB_ENV[t,]))/(1-(Prior_infeB_ENV[t,])*SSe_ENV[t, sample(1:Iter,Iter)])
  
  }
  
  # browser()
  # quantile(Post_freeA[t,],probs=c(0.025,0.5,0.975))
  # quantile(Post_freeA_AFP[t,],probs=c(0.025,0.5,0.975))
  # quantile(Post_freeA_ENV[t,],probs=c(0.025,0.5,0.975))
  
  if (intro){
    # Summarise over iterations, excluding prior column
    free_AFP_sims <- as.data.frame(t(apply(Post_freeB_AFP[-1,],1,quantile,probs=c(0.025,0.5,0.975))))
    names(free_AFP_sims) <- c("lwr","med","upr")
    # print(((0:(max_time-1)/12))+start)
    free_AFP_sims <- mutate(free_AFP_sims,
                            time = ((0:(max_time-1)/12))+start,
                            projection = !is.na(M),
                            future_flag = (time > date_to_period(last_time)))
    
    free_ENV_sims <- as.data.frame(t(apply(Post_freeB_ENV[-1,],1,quantile,probs=c(0.025,0.5,0.975))))
    names(free_ENV_sims) <- c("lwr","med","upr")
    # print(((0:(max_time-1)/12))+start)
    free_ENV_sims <- mutate(free_ENV_sims,
                            time = ((0:(max_time-1)/12))+start,
                            projection = !is.na(M),
                            future_flag = (time > date_to_period(last_time)))
    
    free_AE_sims <-  as.data.frame(t(apply(Post_freeB[-1,],1,quantile,probs=c(0.025,0.5,0.975))))
    names(free_AE_sims) <- c("lwr","med","upr")
    free_AE_sims <- mutate(free_AE_sims,
                           time = ((0:(max_time-1)/12))+start,
                           projection = !is.na(M),
                           future_flag = (time > date_to_period(last_time)))
  }else if(!intro){
    # Summarise over iterations, excluding prior column
    free_AFP_sims <- as.data.frame(t(apply(Post_freeA_AFP[-1,],1,quantile,probs=c(0.025,0.5,0.975))))
    names(free_AFP_sims) <- c("lwr","med","upr")
    # print(((0:(max_time-1)/12))+start)
    free_AFP_sims <- mutate(free_AFP_sims,
                            time = ((0:(max_time-1)/12))+start,
                            projection = !is.na(M),
                            future_flag = (time > date_to_period(last_time)))
    
    free_ENV_sims <- as.data.frame(t(apply(Post_freeA_ENV[-1,],1,quantile,probs=c(0.025,0.5,0.975))))
    names(free_ENV_sims) <- c("lwr","med","upr")
    free_ENV_sims <- mutate(free_ENV_sims,
                            time = ((0:(max_time-1)/12))+start,
                            projection = !is.na(M),
                            future_flag = (time > date_to_period(last_time)))
    
    free_AE_sims <-  as.data.frame(t(apply(Post_freeA[-1,],1,quantile,probs=c(0.025,0.5,0.975))))
    names(free_AE_sims) <- c("lwr","med","upr")
    free_AE_sims <- mutate(free_AE_sims,
                           time = ((0:(max_time-1)/12))+start,
                           projection = !is.na(M),
                           future_flag = (time > date_to_period(last_time)))
  }

  # ---------------------------------------------------------------------------- #
  # Make plots
  
  time_AFP_FFI95 <- min(free_AFP_sims$time[free_AFP_sims$med >= elim_criterion])
  future_flag <- min(free_AFP_sims$time[free_AFP_sims$future])
  p1 <- ggplot(free_AFP_sims,aes(x=time,y=med)) +
    geom_rect(aes(xmin = time_AFP_FFI95,
                  xmax = round(ifelse(is.infinite(time_AFP_FFI95), Inf, (max_time/12)+start)),
                  ymin = 0,
                  ymax = 1),
              fill = "lightblue2",
              alpha = 0.01) +
    geom_rect(aes(xmin = future_flag - 0.03,
                  xmax = round(ifelse(is.infinite(future_flag), Inf, (max_time/12)+start)),
                  ymin = 0,
                  ymax = 1),
              fill = "grey",
              alpha = 0.01) +
    geom_text(aes(x = future_flag + 0.07,y = 0.1), label = "Projection",col="grey", cex = 3) +
    geom_errorbar(aes(ymin=lwr, ymax=upr), width=.01,
                  position=position_dodge(.9),col="grey50") +
    geom_point() +
    geom_hline(aes(yintercept=0.90),linetype="dashed",col="grey50") +
    geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
    geom_text(aes(x = min(time),y = 0.88), label = "0.90",col="grey50", cex = 2) +
    geom_text(aes(x = min(time),y = 0.97), label = "0.99",col="red", cex = 2) +
    labs(x = "Month", y = "Pr(Infection free)") +
    ylim(0,1) +
    ggtitle('AFP Surveillance only')
  
  time_AE_FFI95 <- min(free_AE_sims$time[free_AE_sims$med >= elim_criterion])
  p2 <- ggplot(free_AE_sims,aes(x=time,y=med)) + 
    geom_rect(aes(xmin = time_AE_FFI95,
                  xmax = round(ifelse(is.infinite(time_AE_FFI95), Inf, (max_time/12)+start)),
                  ymin = 0,
                  ymax = 1),
              fill = "lightblue2",
              alpha = 0.01) +
    geom_rect(aes(xmin = future_flag - 0.03,
                  xmax = round(ifelse(is.infinite(future_flag), Inf, (max_time/12)+start)),
                  ymin = 0,
                  ymax = 1),
              fill = "grey",
              alpha = 0.01) +
    geom_text(aes(x = future_flag + 0.07,y = 0.1), label = "Projection",col="grey", cex = 3) +
    geom_errorbar(aes(ymin=lwr, ymax=upr), width=.01,
                  position=position_dodge(.9),col="grey50") +
    geom_point() +
    geom_hline(aes(yintercept=0.95),linetype="dashed",col="grey50") +
    geom_hline(aes(yintercept=0.99),linetype="dashed",col="red") +
    geom_text(aes(x = min(time),y = 0.93), label = "0.90",col="grey50", cex = 2) +
    geom_text(aes(x = min(time),y = 0.97), label = "0.99",col="red", cex = 2) +
    labs(x = "Month", y = "Pr(Infection free)",
         caption = paste("Shaded area in which median FFI probability exceeds",
                         elim_criterion)
    ) +
    ylim(0,1) +
    ggtitle('AFP & ENV Surveillance')
  
  p1/p2
  
  # Define outputs
  out_ffi <- list(SSe_AFP = SSe_AFP,
                  SSe_ENV = SSe_ENV,
                  SSe_AE = SSe_AE,
                  free_AFP_sims = free_AFP_sims,
                  free_ENV_sims = free_ENV_sims,
                  free_AE_sims = free_AE_sims,
                  p1 = p1,
                  p2 = p2)
  
}
