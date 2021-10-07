#### function for simulating survival times with Weibull baseline hazard ####
# eta: vector containing the log hazard (linear predictor) for each subject
# lambda: scale parameter for the baseline hazard
# p: shape paramter for the baseline hazard
# gen_Ct: function for generating censoring times, default is exponential with mean 5 truncated at t=1
gen_St <- function(eta, lambda, p, gen_Ct = function(N) pmin(1,rexp(N, 1/5)) ){
  N <- length(eta)
  U <- runif(N, 0, 1)
  St = (-log(U)/(lambda*exp(eta)))^(1/p)
  Ct = gen_Ct(N)
  
  data.frame("event" = as.numeric(St<=Ct), "time" = pmin(St,Ct), "eta" = eta)
}

#### Integration functions ####
### Note that we need two different functions since there is a double integral in the specificity formula
## function for integrating incident sensitivity
my_fn_sens <- function(x, t, lambda=1, p=1, mn_eta=0, sigma_eta=2){
  lambda*p*t^(p-1)*exp(x[1])*exp(-lambda*t^p*exp(x[1]))*dnorm(x[1], mean=mn_eta, sd=sigma_eta)
}
## function for integrating dynamic specificity
my_fn_spec <- function(x, lambda=1, p=1, mn_eta=0, sigma_eta=2){
  lambda*p*x[2]^(p-1)*exp(x[1])*exp(-lambda*x[2]^p*exp(x[1]))*dnorm(x[1], mean=mn_eta, sd=sigma_eta)
}


#### Calculate AUC ####
## calculate AUC by integrating the ROC curve 
trap_integrate_ROC <- function(eta, sens, spec){
  if(!Inf %in% eta){
    sens <- c(sens, 0)
    spec <- c(spec, 1)
  }
  if(!c(-Inf) %in% eta){
    sens <- c(1, sens)
    spec <- c(0, spec)
  }
  x <- rev(1-spec)
  y <- rev(sens)
  n <- length(y)
  sum( (x[2:n] - x[1:(n-1)] ) * (y[1:(n-1)] + y[2:n]) )/2
}


##### Empirical estimator #####
ID_AUC <- function(marker, Stime, status, predict.time, entry = NULL, ...){
  if (length(entry) == 0) {
    entry = rep(0, NROW(Stime))
  }
  at.risk <- ((Stime >= predict.time) & (entry <= predict.time))
  eta     <- marker[at.risk]
  status2 <- status
  status2[Stime > predict.time] <- 0
  status2 <- status2[at.risk]
  
  C_num <- 0
  n_case_t    <- sum(status2)
  n_control_t <- sum(1-status2)
  inx_ti <- which(status2 == 1)
  inx_tj <- which(status2 == 0)
  for(id in inx_ti){
    C_num <- C_num + sum(eta[id] > eta[inx_tj]) + 0.5*sum(eta[id] == eta[inx_tj])
  }
  out <- C_num/(n_case_t*n_control_t)
  
  out     
}


#### caclulate empirical concordance #####

calc_c <- function(marker, Stime, status){
  utimes <- sort(unique(Stime[status==1]))
  num <- denom <- 0
  for(ut in seq_along(utimes)){
    ## current time
    ti    <- utimes[ut]
    ## subjects who experienced an event at current time
    inx_i <- which(Stime == ti & status==1)
    ## subjects with observed times beyond event current time
    inx_j <- which(Stime > ti)
    ## number of "cases" and "controls" at time current time
    n_case_t    <- length(inx_i)
    n_control_t <- length(inx_j)
    for(i in seq_along(inx_i)){
      num   <- num + sum( (marker[inx_j] > marker[inx_i[i]] ) ) + 0.5*sum(marker[inx_j] == marker[inx_i[i]])
    }
    denom <- denom + n_case_t*n_control_t
  }
  1-num/denom
}

#### calculate concordance weighted by survival probablity ####
intAUC <- function(AUC, utimes, St, method="HZ", smoothAUC=FALSE, n_event=NULL, Ct=NULL, k=30,...){
  ut <- utimes
  # estimate survival probablity
  if(method == "HZ"){
      ft <- rep(NA, length(St))
      ft[1] <- (1 - St[1])/(ut[1])
      for(j in 2:length(St)){
        ft[j] <- (St[j - 1] - St[j])/(ut[j]-ut[j-1])
      }
    }
  if(method == "smS"){
      fit_S <- scam(St ~ s(ut, bs="mpd",k=k),
                    data=data.frame(ut=ut, St=St))
      df_pred <- data.frame(ut = rep(ut, each=2) + rep(c(-0.00001,0.00001), length(ut)))
      St_pred <- predict(fit_S, newdata=df_pred, type='response')
      ft     <- -diff(St_pred)[seq(1,2*length(ut),by=2)]                    
    }
    mIndex <- length(ut)
    wt <- 2 * ft * St
    # use trapezoidal rule to approximate integral
    iAUC <- ut[1] * wt[1] * AUC[1] / 2 # ft[0] = 1 - St[0] = 1 - 1 = 0, so wt[0] = 0
    W <- ut[1] * wt[1] / 2
    for(j in 2:length(St)){
      iAUC <- iAUC + (wt[j] * AUC[j] + wt[j - 1] * AUC[j - 1]) * (ut[j] - ut[j - 1]) / 2
      W <- W + (wt[j] + wt[j - 1]) * (ut[j] - ut[j - 1]) / 2
    }
    iAUC <- iAUC / W
  
  return(iAUC)
}

##### Cauculate true concordance #####

### marginal joint survival function

true_st <- function(t, lambda=2, p=2, eta,sigma_eta){
  St <- exp(-(lambda*exp(eta)*t)^p)*dnorm(eta, 0, sigma_eta)
  return(St)
}

### marginal jpint density function

true_ft <- function(t, lambda=2, p=2, eta,sigma_eta){
  ft <- p*lambda^p*exp(eta*p)*t^(p-1)*exp(-(lambda*exp(eta)*t)^p)*dnorm(eta, 0, sigma_eta)
  return(ft)
}


### distributio of eta: N(0, beta%*%sigma%*%beta)