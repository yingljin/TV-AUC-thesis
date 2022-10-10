
# This script includes functions that are used specifically for data application
# The functions are different from simulation
# mostly because gam model is used instead of Cox PH

##### GH for gam model ##### 

## calculate concordance under Gonen and heller 
## beta_hat: coefficient estimate from cox gam model
## X: design matrix
calc_c_gh <- function(beta_hat, X){
  N <- nrow(X)
  C <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      eta_ij <- (X[i,,drop=F] - X[j,,drop=F]) %*% beta_hat
      if(eta_ij < 0){
        C <- C + 1/(1 + exp(eta_ij))
      }
      else if(eta_ij>0){
        C <- C + 1/(1 + exp(-eta_ij))
      }
     
   }
  }
  2/(N*(N-1))*C
}


#### Training estimates #### 

# calcualte AUC from GAM model on training set

## fit: GAM-cox model
## data: dataframe with columns time, event
## t: unique event time in training data
## nt: number of unique event time


train_auc <- function(fit = fit_k, data = df_train_k,
                      t = t_uni_train, nt = nt_uni_train){
  # set up results container
  auc_mat <- matrix(NA, nrow = nt, ncol = 4)
  colnames(auc_mat) <- c("time", "HZ", "NP", "SNP")
  auc_mat[, "time"] <- t
  
  # estimated risk score
  data$eta <- fit$linear.predictors
  
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
  auc_mat[, "SNP"] <- predict(sm_fit)
  
  return(auc_mat)
}


#### Testing estimates ####

## eta: linear predictors/risk score
## data: dataframe with columns time, event
## t: unique event time in testing data
## nt: number of unique event time

test_auc <- function(eta, data = data_test, t = t_uni_test, nt = nt_uni_test){
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
  sm_fit <- gam(NP ~ s(time, bs = "cr", k = 20), method = "REML", 
                data = as.data.frame(auc_mat))
  auc_mat[, "SNP"] <- predict(sm_fit)
  
  return(auc_mat)
}


#### Concrodance from GAM-Cox model ####

# calculate truncated concordance weighted by estimated survival probability and density
# AUC: AUC estimates from previous functions
# utimes: event time
# St: survival function values
# method: "HZ" estimates survival density by slope, while "smS" uses scam smoothing


intAUC_appl <- function(AUC, utimes, St, method="HZ"){
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
    fit_S <- scam(St ~ s(ut, bs="mpd", k = 20),
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
