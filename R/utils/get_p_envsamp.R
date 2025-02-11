get_p_envsamp <- function(freq_est, Iter = 1000){
  
  samp_vals <- read_csv(here("inputs","ES_envsamp_vals.csv"))
  draws <- matrix(ncol = length(freq_est), nrow = Iter)
  
  for (i in seq_along(freq_est)){
    
    if (is.na(freq_est[i]) | freq_est[i] == 0){
      # print(paste0("Missing frequency estimate in ",i))
      draws[,i] <- 0
    }else{
      # Parameters cannot equal exactly 1
      tmp <- betaExpert(best=min(samp_vals$EstMean[samp_vals$SampleFreq == max(round(freq_est[i]),1)],0.999),
                        lower=min(samp_vals$EstLwr[samp_vals$SampleFreq == max(round(freq_est[i]),1)],0.999),
                        upper=min(samp_vals$EstUpr[samp_vals$SampleFreq == max(round(freq_est[i]),1)],0.999),
                        p=0.95)
      
      draws[,i] <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta) 
    }

  }

  return(draws)
  
}
