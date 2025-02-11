period_to_date = function(x,format='%Y-%b'){
  yr = floor(x)
  mo = round(12*(x-yr),0)+1
  d = '1'
  format(ymd(paste(yr,mo,d,sep='-')),format)
}

date_to_period = function(x){
  year(x) + (month(x)-1)/12
}

#### Create radiation matrix ##########

radiation_matrix_construction = function(x,pop,threshold = 1e-6){
  deg2rad <- function(deg) return(deg*pi/180)
  e_radius <- 6371 # Earth mean radius [km]
  
  D = x %>% as.matrix %>% deg2rad
  D <- acos(outer(sin(D[,2]),sin(D[,2])) + outer(cos(D[,2]),cos(D[,2]))* cos(outer(D[,1],D[,1],FUN='-'))) *e_radius
  diag(D) = 0
  R = matrix(0,nrow = nrow(D),ncol=ncol(D))
  
  for(i in 1:nrow(D)){
    ii = order(D[i,])
    p =pop[ii]
    s = cumsum(p)
    R[i,ii] = (p[1]*p)/((s-p)*(s))
    R[i,i] <- 1
  }
  
  R = Matrix(R)
  R[is.na(R)] = 0
  R[is.infinite(R)] = 0
  
  #make a sparse matrix to cut down on memory/time
  R[R< threshold] = 0
  R = as(R, "sparseMatrix")
  
  R
}

#### Get pr sampling during shedding ##########

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


#### Solve for adjusted risks ##########

ARfunct2 <- function(data) {
  # add in some check
  if(round(sum(data$PrP),10)!=1){
    warning("check sum of PrP")
  }
  if(sum(data$RR==1)<1){
    warning("check that 1 RR == 1.00")
  }
  RR <- data$RR          # estimated risk
  PrP <- data$PrP
  n <- length(RR)                          # number of equations
  # which RR==1
  small <- which(RR==1)
  AR <- rep(NA,n)
  # What is AR_small? ie. eqn 1
  AR_tmp <- rep(NA,n-length(small))
  tmp1 <- sum(PrP[small])
  tmp2 <- RR[-small]*PrP[-small]
  AR1 <- 1/(tmp1 + sum(tmp2))
  # eqn 2 - use AR1 in eqns
  vals <- c(1:n)
  vals2 <- vals[-small]
  AR[small] <- AR1  # those with RR=1
  for(i in vals2){
    AR[i] <- RR[i]*AR1  # those left
  }
  #browser()
  return(AR)
} 