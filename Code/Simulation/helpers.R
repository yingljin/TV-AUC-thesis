
# This script includes functions written for data generation
# and calculating the estimators 


#### data generation  ####
# function for simulating survival times with Weibull baseline hazard
# eta: vector containing the linear predictor for each subject
# lambda: scale parameter for the baseline hazard
# p: shape paramter for the baseline hazard
# gen_Ct: function for generating censoring times
gen_St <- function(eta, lambda=2, p=2, gen_Ct){
  N <- length(eta)
  U <- runif(N, 0, 1)
  St = (-log(U)/(lambda*exp(eta)))^(1/p)
  Ct = gen_Ct(N)
  
  data.frame("event" = as.numeric(St<=Ct), "time" = pmin(St,Ct), "eta" = eta)
}


##### Non-parametric estimator of AUC #####
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


#### Harrell's C-index #####
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

#### Truncated Concordance ####
#  integral of AUC weighted by smoothed estimated survival probability and density

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

##### Gonen & Heller #####
coxphCPE_eta <- function (betahat,Xmat){
  # if (class(phfit) != "coxph") 
  #   stop("phfit shoud be coxph class object")
  n <- nrow(Xmat)
  # betahat <- phfit$coefficients
  p <- length(betahat)
  # vbetahat <- phfit$var
  xbeta <- Xmat%*%betahat
  bw <- 0.5 * sd(xbeta) * (n^(-1/3))
  zzz <- .Fortran("cpesub", as.integer(n), as.integer(p), as.double(Xmat), 
                  as.double(xbeta), as.double(bw), CPE = double(1), CPEsmooth = double(1), 
                  varDeriv = double(p), uRowSum = double(n), uSSQ = double(1), 
                  PACKAGE = "clinfun")
  CPE <- 2 * zzz$CPE/(n * (n - 1))
  # CPEsmooth <- 2 * zzz$CPEsmooth/(n * (n - 1))
  # varTerm1 <- 4 * (sum((zzz$uRowSum + rep(0.5, n) - n * CPEsmooth)^2) -
  #                    (zzz$uSSQ + n/4 - n * CPEsmooth - n * (n - 2) * CPEsmooth^2))/(n *
  #                                                                                     (n - 1))^2
  # varDeriv <- 2 * zzz$varDeriv/(n * (n - 1))
  # varTerm2 <- t(varDeriv) %*% vbetahat %*% varDeriv
  # varCPE <- varTerm1 + varTerm2
  # out <- c(CPE, CPEsmooth, sqrt(varCPE))
  # names(out) <- c("CPE", "smooth.CPE", "se.CPE")
  return(CPE)
}


