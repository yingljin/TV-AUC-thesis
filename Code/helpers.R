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
