# This script includes functions for aggregating results 
# for one simulated dataset

#### In-sample estimates  #### 

# AUC
## fit is the cox model
## t: unique event time in training data
## eta: risk score corresponding to t 

train_auc <- function(fit, data = data_train,
                      t = t_uni_train, nt = nt_uni_train){
  # set up results container
  auc_mat <- matrix(NA, nrow = nt, ncol = 4)
  colnames(auc_mat) <- c("time", "HZ", "NP", "SNP")
  auc_mat[, "time"] <- t
  
  # estimated risk score
  data$eta <- fit$linear.predictors
  
  # cauclate in-sample AUC
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

# concordance 
train_concord <- function(data = data_train, KM_est = KM_est_train, 
                          auc_mat, fit){
  
  auc_sort <-arrange(data.frame(auc_mat), time)
  
  # result container
  c_mat <- matrix(NA, nrow = 1, ncol = 5)
  colnames(c_mat)<- c("GH", "HZ",
                      "Harrell's", "NP", "SNP")
  # concordance
  c_mat[, "GH"] <- coxphCPE(fit)["CPE"]
  c_mat[, "HZ"] <- intAUC(auc_sort$HZ, auc_sort$time, 
                          KM_est, method = "HZ")
  c_mat[, "Harrell's"] <- calc_c(fit$linear.predictors, data$time, data$event)
  c_mat[, "NP"] <- intAUC(auc_sort$NP, auc_sort$time, 
                          KM_est, method = "smS")
  c_mat[, "SNP"] <- intAUC(auc_sort$SNP, auc_sort$time, 
                           KM_est, method = "smS")
  
  return(c_mat)
}


#### Out-of-sample estimates ####

# eta is a vector biomaker/linear predictor

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
  sm_fit <- gam(NP ~ s(time, k = 30, bs = "cr"), method = "REML", 
                data = as.data.frame(auc_mat))
  auc_mat[, "SNP"] <- predict(sm_fit)
  
  return(auc_mat)
}


# test concordance
test_concord <- function(data = data_test, KM_est = KM_est_test, 
                         auc_mat, fit, eta, X_mat){
  
  auc_sort <-arrange(data.frame(auc_mat), time)
  
  # result container
  c_mat <- matrix(NA, nrow = 1, ncol = 5)
  colnames(c_mat)<- c("GH", "HZ",
                      "Harrell's", "NP", "SNP")
  # concordance
  c_mat[, "GH"] <- coxphCPE_eta(fit, X_mat)['CPE']
  c_mat[, "HZ"] <- intAUC(auc_sort$HZ, auc_sort$time, 
                          KM_est, method = "HZ")
  c_mat[, "Harrell's"] <- calc_c(eta, data$time, data$event)
  c_mat[, "NP"] <- intAUC(auc_sort$NP, auc_sort$time, 
                          KM_est, method = "smS")
  c_mat[, "SNP"] <- intAUC(auc_sort$SNP, auc_sort$time, 
                           KM_est, method = "smS")
  
  return(c_mat)
}
