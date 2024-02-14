
# This script includes functions written for data generation
# and calculating the estimators 


#### data generation  ####
# function for simulating survival times with Weibull baseline hazard
# eta: vector containing the linear predictor for each subject
# lambda: scale parameter for the baseline hazard
# p: shape parameter for the baseline hazard
# gen_Ct: function for generating censoring times

gen_St <- function(eta, lambda=2, p=2, gen_Ct){
  N <- length(eta)
  U <- runif(N, 0, 1)
  St = (-log(U)/(lambda*exp(eta)))^(1/p)
  Ct = gen_Ct(N)
  
  data.frame("event" = as.numeric(St<=Ct), "time" = pmin(St,Ct), "eta" = eta)
}


#### Non-parametric AUC ####
# makers: estimated risk score/linear predictors
# Stime: even times
# status: even status
# predict.time: time point at which to calculate AUC for 

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


#### Concrodance by integrating AUC ####

# calculate truncated concordance weighted by estimated survival probability and density
# AUC: AUC estimates from previous functions
# utimes: event time
# St: survival function values
# method: "HZ" estimates survival density by slope, while "smS" uses scam smoothing

intAUC <- function(AUC, utimes, St, method="HZ", k=10,...){
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

#### Harrell's C-index #####
# makers: estimated risk score/linear predictors
# Stime: even times
# status: even status

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
##### Gonen & Heller #####
# betahat: estimates of coefficient
# Xmat: covariates 

coxphCPE_eta <- function (betahat,Xmat){
  n <- nrow(Xmat)
  p <- length(betahat)
  xbeta <- Xmat%*%betahat
  bw <- 0.5 * sd(xbeta) * (n^(-1/3))
  zzz <- .Fortran("cpesub", as.integer(n), as.integer(p), as.double(Xmat), 
                  as.double(xbeta), as.double(bw), CPE = double(1), CPEsmooth = double(1), 
                  varDeriv = double(p), uRowSum = double(n), uSSQ = double(1), 
                  PACKAGE = "clinfun")
  CPE <- 2 * zzz$CPE/(n * (n - 1))
  return(CPE)
}

#### Calculate and format AUC ####

# eta is a vector biomaker/linear predictor
## t: unique event time in training data
## nt: number of unique events time points
## data: data set with columns of time, event


tv_auc <- function(eta, data, t, nt){
  # set up results container
  auc_mat <- matrix(NA, nrow = nt, ncol = 4)
  colnames(auc_mat) <- c("time", "HZ", "NP", "SNP")
  auc_mat[, "time"] <- t
  
  data$eta <- eta
  
  # calculate in-sample AUC
  for(i in seq_along(t)){
    t_i <- t[i]
    # HZ
    auc_mat[i, "HZ"] <- CoxWeights(marker = data$eta, 
                                   Stime = data$time,
                                   status = data$event,
                                   predict.time = t_i, entry = NULL)$AUC
    
    # non-parametric estimator
    auc_mat[i, "NP"] <- ID_AUC(marker = data$eta, 
                               Stime = data$time, 
                               status = data$event,
                               predict.time = t_i,  entry = NULL)
    
  }
  # smoothed empirical AUC
  sm_fit <- gam(NP ~ s(time, k = 30, bs = "cr"), method = "REML", 
                data = as.data.frame(auc_mat))
  auc_mat[, "SNP"] <- predict(sm_fit, newdata = data.frame(time = auc_mat[, "time"]))
  
  return(auc_mat)
}


#### Calculate and format AUC ####

## data: data set with columns of time, event
## KM_est: estimated survival 
## auc_mat: time-varying auc 
## coef: estimated coefficients
## X_mat: covariantes

concord <- function(data, KM_est, auc_mat, coef, X_mat){
  
  auc_sort <-arrange(data.frame(auc_mat), time)
  eta <- X_mat %*% coef
  
  # result container
  c_mat <- matrix(NA, nrow = 1, ncol = 5)
  colnames(c_mat)<- c("GH", "HZ",
                      "Harrell's", "NP", "SNP")
  # concordance
  c_mat[, "GH"] <- coxphCPE_eta(betahat = coef, X_mat)
  c_mat[, "HZ"] <- intAUC(auc_sort$HZ, auc_sort$time, 
                          KM_est, method = "HZ")
  c_mat[, "Harrell's"] <- calc_c(eta, data$time, data$event)
  c_mat[, "NP"] <- intAUC(auc_sort$NP, auc_sort$time, 
                          KM_est, method = "smS")
  c_mat[, "SNP"] <- intAUC(auc_sort$SNP, auc_sort$time, 
                           KM_est, method = "smS")
  
  return(c_mat)
}

