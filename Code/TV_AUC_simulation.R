rm(list = ls())
library("survival")
library("ggplot2")
library("cubature")
library("tidyverse")
library(here)
library(risksetROC)
library(mgcv)
library(scam)
library("clinfun")
# set the seed
set.seed(102131)


# helper functions
source(here("Code/helpers.R"))

#### Simulation set up ####

M <- 1000 # number of simulated data sets
N <- 500 # number of subjects
Beta <- c(1,-1,0.25)

# training sample
N_obs <- N*0.5

# delta_t <- 0.05

# result container
auc_lst_train <- list()
auc_lst_test <- list()
c_lst_train <- list()
c_lst_test <- list()
data_lst <- list()

##### simulation ####
pb <- txtProgressBar(min=0, max=M,style=3)
for(iter in 1:M){
  # generate data
  X  <- matrix(rnorm(N*3), ncol=3, nrow=N)
  ## add noise signal: 20 or 100, N(0, 1)
  Z <- matrix(rnorm(N*100), ncol = 100, nrow = N)
  eta <- X %*% Beta
  data <- gen_St(eta=eta,lambda=2, p=2, gen_Ct = function(N){sample(c(0.5, 1), size = N, replace = T)})
  data$X <- I(X)
  data$Z <- I(Z)
  data_lst[[iter]] <- data
  ## separate out the training and test datasets
  data_train <- data[1:N_obs,]
  data_test  <- data[-c(1:N_obs),]
  
  ## fit our 3 models with different number of noise predictors (0, 20, 100)
  fit_1 <- coxph(Surv(time, event) ~ X, data=data_train)
  fit_2 <- coxph(Surv(time, event) ~ X + Z[,1:20], data=data_train)
  fit_3 <- coxph(Surv(time, event) ~ X + Z, data=data_train)
  
  # get estimated \hat\eta_ik for k = 1,2,3 models, unique time and biomarker
  ## for the training data
  data_train$eta1 <- fit_1$linear.predictors
  data_train$eta2 <- fit_2$linear.predictors
  data_train$eta3 <- fit_3$linear.predictors
  t_vec_train <- unique(data_train$time[data_train$event==1])
  eta_vec_train <- unique(data_train$eta[data_train$event==1])
  nt_pred_train <- length(t_vec_train)
  neta_pred_train <- length(eta_vec_train)
  
  ## for the test data
  data_test$eta1 <- as.vector(data_test$X %*% coef(fit_1))
  data_test$eta2 <- as.vector(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2))
  data_test$eta3 <- as.vector(cbind(data_test$X,data_test$Z) %*% coef(fit_3))
  t_vec_test <- unique(data_test$time[data_test$event==1])
  eta_vec_test <- unique(data_test$eta[data_test$event==1])
  nt_pred_test <- length(t_vec_test)
  neta_pred_test <- length(eta_vec_test)
  
  # set up results container
  auc_mat_train <- matrix(NA, nrow = nt_pred_train, ncol = 3*3+1)
  auc_mat_test <- matrix(NA, nrow = nt_pred_test, ncol = 3*3+1)
  colnames(auc_mat_train)<-colnames(auc_mat_test) <- c("time", "HZ_eta1", "HZ_eta2", "HZ_eta3", 
                         "empirical_eta1", "empirical_eta2", "empirical_eta3",
                         "sm_empirical_eta1", "sm_empirical_eta2", "sm_empirical_eta3")
  auc_mat_train[, "time"] <- t_vec_train    
  auc_mat_test[, "time"] <- t_vec_test   
  
  c_mat_train <- matrix(NA, nrow = 3, ncol = 8)
  c_mat_test <- matrix(NA, nrow = 3, ncol = 8)
  colnames(c_mat_train)<-colnames(c_mat_test) <- c("empirical", "GH", "HZ_HZ", "empirical_HZ", "sm_empirical_HZ",
                       "HZ_SmS", "empirical_SmS", "sm_empirical_SmS")
  rownames(c_mat_train)<-rownames(c_mat_test)  <- c("eta1", "eta2", "eta3")
  
  # calculate in sample HZ and empirical time_varying AUC
  for(i in seq_along(t_vec_train)){
    t_i <- t_vec_train[i]
    # HZ
    auc_mat_train[i, "HZ_eta1"] <- CoxWeights(marker = data_train$eta1, 
                                              Stime = data_train$time, status = data_train$event,
                                      predict.time = t_i, entry = NULL)$AUC
    auc_mat_train[i, "HZ_eta2"] <- CoxWeights(marker = data_train$eta2, 
                                              Stime = data_train$time, status = data_train$event,
                                              predict.time = t_i, entry = NULL)$AUC
    auc_mat_train[i, "HZ_eta3"] <- CoxWeights(marker = data_train$eta3, 
                                              Stime = data_train$time, status = data_train$event,
                                              predict.time = t_i, entry = NULL)$AUC
    # trying to get non-parametric estimator
    auc_mat_train[i, "empirical_eta1"] <- ID_AUC(marker = data_train$eta1, 
                                           Stime = data_train$time, status = data_train$event,
                                  predict.time = t_i,  entry = NULL)
    auc_mat_train[i, "empirical_eta2"] <- ID_AUC(marker = data_train$eta2, 
                                                 Stime = data_train$time, status = data_train$event,
                                                 predict.time = t_i,  entry = NULL)
    auc_mat_train[i, "empirical_eta3"] <- ID_AUC(marker = data_train$eta3, 
                                                 Stime = data_train$time, status = data_train$event,
                                                 predict.time = t_i,  entry = NULL)
  }
  # smoothed empirical AUC
  sm_fit1_train <- gam(empirical_eta1 ~ s(time, k = 30, bs = "cr"), method = "REML", 
                 data = as.data.frame(auc_mat_train))
  sm_fit2_train <- gam(empirical_eta2 ~ s(time, k = 30, bs = "cr"), method = "REML", 
                       data = as.data.frame(auc_mat_train))
  sm_fit3_train <- gam(empirical_eta3 ~ s(time, k = 30, bs = "cr"), method = "REML", 
                       data = as.data.frame(auc_mat_train))
  
  auc_mat_train[, "sm_empirical_eta1"] <- predict(sm_fit1_train)
  auc_mat_train[, "sm_empirical_eta2"] <- predict(sm_fit2_train)
  auc_mat_train[, "sm_empirical_eta3"] <- predict(sm_fit3_train)
  
  # calculate out-of-sample HZ and empirical time_varying AUC
  for(i in seq_along(t_vec_test)){
    t_i <- t_vec_test[i]
    # HZ
    auc_mat_test[i, "HZ_eta1"] <- CoxWeights(marker = data_test$eta1, 
                                              Stime = data_test$time, status = data_test$event,
                                              predict.time = t_i, entry = NULL)$AUC
    auc_mat_test[i, "HZ_eta2"] <- CoxWeights(marker = data_test$eta2, 
                                              Stime = data_test$time, status = data_test$event,
                                              predict.time = t_i, entry = NULL)$AUC
    auc_mat_test[i, "HZ_eta3"] <- CoxWeights(marker = data_test$eta3, 
                                              Stime = data_test$time, status = data_test$event,
                                              predict.time = t_i, entry = NULL)$AUC
    # trying to get non-parametric estimator
    auc_mat_test[i, "empirical_eta1"] <- ID_AUC(marker = data_test$eta1, 
                                                 Stime = data_test$time, status = data_test$event,
                                                 predict.time = t_i,  entry = NULL)
    auc_mat_test[i, "empirical_eta2"] <- ID_AUC(marker = data_test$eta2, 
                                                 Stime = data_test$time, status = data_test$event,
                                                 predict.time = t_i,  entry = NULL)
    auc_mat_test[i, "empirical_eta3"] <- ID_AUC(marker = data_test$eta3, 
                                                 Stime = data_test$time, status = data_test$event,
                                                 predict.time = t_i,  entry = NULL)
  }
  # smoothed empirical AUC
  sm_fit1_test <- gam(empirical_eta1 ~ s(time, k = 30, bs = "cr"), method = "REML", 
                       data = as.data.frame(auc_mat_test))
  sm_fit2_test <- gam(empirical_eta2 ~ s(time, k = 30, bs = "cr"), method = "REML", 
                       data = as.data.frame(auc_mat_test))
  sm_fit3_test <- gam(empirical_eta3 ~ s(time, k = 30, bs = "cr"), method = "REML", 
                       data = as.data.frame(auc_mat_test))
  
  auc_mat_test[, "sm_empirical_eta1"] <- predict(sm_fit1_test)
  auc_mat_test[, "sm_empirical_eta2"] <- predict(sm_fit2_test)
  auc_mat_test[, "sm_empirical_eta3"] <- predict(sm_fit3_test)
  
  # concordance
  ## First, estimate survival probability from Kaplan Meier curve 
  KM_fit_train <- survfit(Surv(time,event)~1, timefix=FALSE,data=data_train)
  KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]
  KM_fit_test <- survfit(Surv(time,event)~1, timefix=FALSE,data=data_test)
  KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]
  
  
  auc_sort_train <-arrange(data.frame(auc_mat_train), time)
  auc_sort_test <-arrange(data.frame(auc_mat_test), time)
  ## in-sample HZ and smS concordance
  c_mat_train["eta1", "HZ_HZ"] <- intAUC(auc_sort_train$HZ_eta1, auc_sort_train$time, 
                                         KM_est_train, method = "HZ")
  c_mat_train["eta2", "HZ_HZ"] <- intAUC(auc_sort_train$HZ_eta2, auc_sort_train$time, 
                                         KM_est_train, method = "HZ")
  c_mat_train["eta3", "HZ_HZ"] <- intAUC(auc_sort_train$HZ_eta3, auc_sort_train$time, 
                                         KM_est_train, method = "HZ")
  c_mat_train["eta1", "empirical_HZ"] <- intAUC(auc_sort_train$empirical_eta1, auc_sort_train$time, 
                                                KM_est_train, method = "HZ")
  c_mat_train["eta2", "empirical_HZ"] <- intAUC(auc_sort_train$empirical_eta2, auc_sort_train$time, 
                                                KM_est_train, method = "HZ")
  c_mat_train["eta3", "empirical_HZ"] <- intAUC(auc_sort_train$empirical_eta3, auc_sort_train$time, 
                                                KM_est_train, method = "HZ")
  c_mat_train["eta1", "sm_empirical_HZ"] <- intAUC(auc_sort_train$sm_empirical_eta1, auc_sort_train$time,
                                                   KM_est_train, method = "HZ")
  c_mat_train["eta2", "sm_empirical_HZ"] <- intAUC(auc_sort_train$sm_empirical_eta2, auc_sort_train$time,
                                                   KM_est_train, method = "HZ")
  c_mat_train["eta3", "sm_empirical_HZ"] <- intAUC(auc_sort_train$sm_empirical_eta3, auc_sort_train$time,
                                                   KM_est_train, method = "HZ")
  c_mat_train["eta1", "HZ_SmS"] <- intAUC(auc_sort_train$HZ_eta1, auc_sort_train$time, 
                                          KM_est_train, method = "smS")
  c_mat_train["eta2", "HZ_SmS"] <- intAUC(auc_sort_train$HZ_eta2, auc_sort_train$time, 
                                          KM_est_train, method = "smS")
  c_mat_train["eta3", "HZ_SmS"] <- intAUC(auc_sort_train$HZ_eta3, auc_sort_train$time, 
                                          KM_est_train, method = "smS")
  c_mat_train["eta1", "empirical_SmS"] <- intAUC(auc_sort_train$empirical_eta1, auc_sort_train$time, 
                                                 KM_est_train, method = "smS")
  c_mat_train["eta2", "empirical_SmS"] <- intAUC(auc_sort_train$empirical_eta2, auc_sort_train$time, 
                                                 KM_est_train, method = "smS")
  c_mat_train["eta3", "empirical_SmS"] <- intAUC(auc_sort_train$empirical_eta3, auc_sort_train$time, 
                                                 KM_est_train, method = "smS")
  c_mat_train["eta1", "sm_empirical_SmS"] <- intAUC(auc_sort_train$sm_empirical_eta1, auc_sort_train$time, 
                                                    KM_est_train, method = "smS")
  c_mat_train["eta2", "sm_empirical_SmS"] <- intAUC(auc_sort_train$sm_empirical_eta2, auc_sort_train$time, 
                                                    KM_est_train, method = "smS")
  c_mat_train["eta3", "sm_empirical_SmS"] <- intAUC(auc_sort_train$sm_empirical_eta3, auc_sort_train$time, 
                                                    KM_est_train, method = "smS")
  ### Harrels C-index
  c_mat_train["eta1", "empirical"] <- calc_c(data_train$eta1, data_train$time, data_train$event)
  c_mat_train["eta2", "empirical"] <- calc_c(data_train$eta2, data_train$time, data_train$event)
  c_mat_train["eta3", "empirical"] <- calc_c(data_train$eta3, data_train$time, data_train$event)
  ### Gonen-Heller
  c_mat_train["eta1", "GH"] <- coxphCPE(fit_1)["CPE"]
  c_mat_train["eta2", "GH"] <- coxphCPE(fit_2)["CPE"]
  c_mat_train["eta3", "GH"] <- coxphCPE(fit_3)["CPE"]
  
  ## out-of-sample HZ and smS concordance
  c_mat_test["eta1", "HZ_HZ"] <- intAUC(auc_sort_test$HZ_eta1, auc_sort_test$time, 
                                         KM_est_test, method = "HZ")
  c_mat_test["eta2", "HZ_HZ"] <- intAUC(auc_sort_test$HZ_eta2, auc_sort_test$time, 
                                         KM_est_test, method = "HZ")
  c_mat_test["eta3", "HZ_HZ"] <- intAUC(auc_sort_test$HZ_eta3, auc_sort_test$time, 
                                         KM_est_test, method = "HZ")
  c_mat_test["eta1", "empirical_HZ"] <- intAUC(auc_sort_test$empirical_eta1, auc_sort_test$time, 
                                                KM_est_test, method = "HZ")
  c_mat_test["eta2", "empirical_HZ"] <- intAUC(auc_sort_test$empirical_eta2, auc_sort_test$time, 
                                                KM_est_test, method = "HZ")
  c_mat_test["eta3", "empirical_HZ"] <- intAUC(auc_sort_test$empirical_eta3, auc_sort_test$time, 
                                                KM_est_test, method = "HZ")
  c_mat_test["eta1", "sm_empirical_HZ"] <- intAUC(auc_sort_test$sm_empirical_eta1, auc_sort_test$time,
                                                   KM_est_test, method = "HZ")
  c_mat_test["eta2", "sm_empirical_HZ"] <- intAUC(auc_sort_test$sm_empirical_eta2, auc_sort_test$time,
                                                   KM_est_test, method = "HZ")
  c_mat_test["eta3", "sm_empirical_HZ"] <- intAUC(auc_sort_test$sm_empirical_eta3, auc_sort_test$time,
                                                   KM_est_test, method = "HZ")
  c_mat_test["eta1", "HZ_SmS"] <- intAUC(auc_sort_test$HZ_eta1, auc_sort_test$time, 
                                          KM_est_test, method = "smS")
  c_mat_test["eta2", "HZ_SmS"] <- intAUC(auc_sort_test$HZ_eta2, auc_sort_test$time, 
                                          KM_est_test, method = "smS")
  c_mat_test["eta3", "HZ_SmS"] <- intAUC(auc_sort_test$HZ_eta3, auc_sort_test$time, 
                                          KM_est_test, method = "smS")
  c_mat_test["eta1", "empirical_SmS"] <- intAUC(auc_sort_test$empirical_eta1, auc_sort_test$time, 
                                                 KM_est_test, method = "smS")
  c_mat_test["eta2", "empirical_SmS"] <- intAUC(auc_sort_test$empirical_eta2, auc_sort_test$time, 
                                                 KM_est_test, method = "smS")
  c_mat_test["eta3", "empirical_SmS"] <- intAUC(auc_sort_test$empirical_eta3, auc_sort_test$time, 
                                                 KM_est_test, method = "smS")
  c_mat_test["eta1", "sm_empirical_SmS"] <- intAUC(auc_sort_test$sm_empirical_eta1, auc_sort_test$time, 
                                                    KM_est_test, method = "smS")
  c_mat_test["eta2", "sm_empirical_SmS"] <- intAUC(auc_sort_test$sm_empirical_eta2, auc_sort_test$time, 
                                                    KM_est_test, method = "smS")
  c_mat_test["eta3", "sm_empirical_SmS"] <- intAUC(auc_sort_test$sm_empirical_eta3, auc_sort_test$time, 
                                                    KM_est_test, method = "smS")
  ### Harrels C-index
  c_mat_test["eta1", "empirical"] <- calc_c(data_test$eta1, data_test$time, data_test$event)
  c_mat_test["eta2", "empirical"] <- calc_c(data_test$eta2, data_test$time, data_test$event)
  c_mat_test["eta3", "empirical"] <- calc_c(data_test$eta3, data_test$time, data_test$event)
  ### Gonen-Heller: need to use the fitted cox to calculate biomarker
  # fit_1_test <- coxph(Surv(time, event) ~ X, data=data_test)
  # fit_2_test <- coxph(Surv(time, event) ~ X + Z[,1:20], data=data_test)
  # fit_3_test <- coxph(Surv(time, event) ~ X + Z, data=data_test)
  c_mat_test["eta1", "GH"] <- coxphCPE_eta(fit_1, data_test$X)['CPE']
  c_mat_test["eta2", "GH"] <- coxphCPE_eta(fit_2, cbind(data_test$X, data_test$Z[, 1:20]))["CPE"]
  c_mat_test["eta3", "GH"] <-  coxphCPE_eta(fit_3, cbind(data_test$X, data_test$Z))["CPE"]
  
  # save to final results
  auc_lst_train[[iter]] <- auc_mat_train
  auc_lst_test[[iter]] <- auc_mat_test
  c_lst_train[[iter]] <- c_mat_train
  c_lst_test[[iter]] <- c_mat_test
  
  print(iter)
  setTxtProgressBar(pb, value=iter)
}

auc_lst_train[[1]] %>% head()
auc_lst_test[[1]] %>% head()
c_lst_train[[1]]
c_lst_test[[1]]

auc_lst_train <- lapply(auc_lst_train, as.data.frame)
auc_lst_test <- lapply(auc_lst_test, as.data.frame)
auc_df_train <- bind_rows(auc_lst_train, .id = "iter")
auc_df_test <- bind_rows(auc_lst_test, .id = "iter")
length(unique((auc_df_train$iter)))

## concordance

c_lst_train <- lapply(c_lst_train, as.data.frame)
c_df_train <- lapply(c_lst_train, rownames_to_column, "signal") 
c_df_train <- bind_rows(c_df_train, .id = "iter")
c_df_train <- c_df_train %>% pivot_longer(3:10, names_to = "estimator", values_to = "concordance")
c_df_train$estimator <- factor(c_df_train$estimator, 
                         levels = c("empirical", "GH", "HZ_HZ", "HZ_SmS",
                                    "empirical_HZ", "empirical_SmS",
                                    "sm_empirical_HZ", "sm_empirical_SmS"))

c_lst_test <- lapply(c_lst_test, as.data.frame)
c_df_test <- lapply(c_lst_test, rownames_to_column, "signal") 
c_df_test <- bind_rows(c_df_test, .id = "iter")
c_df_test <- c_df_test %>% pivot_longer(3:10, names_to = "estimator", values_to = "concordance")
c_df_test$estimator <- factor(c_df_test$estimator, 
                               levels = c("empirical", "GH", "HZ_HZ", "HZ_SmS",
                                          "empirical_HZ", "empirical_SmS",
                                          "sm_empirical_HZ", "sm_empirical_SmS"))



save(auc_df_train, auc_df_test, c_df_train, c_df_test, file = here("outputData/estimated_values.RData"))
save(auc_lst_train, auc_lst_test, file = here("outputData/results_by_iter.RData"))

##### true AUC and concordance #####
# do not run this section
# instead load saved data 
# unless we need different method to calculated true AUC
range(auc_df$time)
tind <- seq(0, 1, length.out = 500)
etaind <- seq(-5, 5, length.out = 100)
true_auc <- rep(NA, length(tind))

summary(auc_df$time)
summary(data$eta)

## in calculating sensitivity and specificity
## restrict boundaries of eta and t to a reasonable value based on distribution
## to avoid numeric issues
## technically, time cannot be greater than 1
pb <- txtProgressBar(min=0, max=length(tind), style=3)
for(i in seq_along(tind)){
  t_i <- tind[i]
  ## containers
  sens <- rep(NA, length(etaind))
  spec <- rep(NA, length(etaind))
  for(j in seq_along(etaind)){
    eta_ij <- etaind[j]
    
    sens[j] <- adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(eta_ij), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral/
      adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(-Inf), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral
    
    spec[j] <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(eta_ij, 500), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral/
      adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(500, 500), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral
  }
  ## integrate to AUC
  true_auc[i] <- trap_integrate_ROC(etaind, sens, spec)
  
  setTxtProgressBar(pb, value=i)
}

# brief look at true auc
plot(tind, true_auc) 

# true concordance
## didn't include the interpolated AUC
## because they are technically estimated
true_auc_sort <- data.frame(time_bin = tind, auc = true_auc)

## shall we use marginal survival function by intergrating out eta?
## marginal S(t)
sig_eta <- sqrt(sum(Beta^2))
true_marg_st <- rep(NA, length(tind))
for(i in seq_along(true_marg_st)){
  true_marg_st[i] <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, tind[i]), upperLimit=c(Inf, 500), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral
}

plot(tind, true_marg_st)

## marginal f(t)
true_marg_ft <- rep(NA, length(tind))
for(i in seq_along(true_marg_ft)){
  true_marg_ft[i] <- adaptIntegrate(my_fn_sens, t=tind[i], lowerLimit=c(-Inf), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral
}

plot(tind, true_marg_ft)


## use trapezoidal rule to approximate integral
## use truncated version for weights
w_vec <- 2 * true_marg_ft * true_marg_st
w_vec <- w_vec/(1-true_marg_st[which.max(tind)]^2)
y_vec <- w_vec*true_auc_sort$auc

plot(tind, w_vec)
plot(tind, y_vec)

nt <- nrow(true_auc_sort)
width <- diff(true_auc_sort$time_bin)
height <- y_vec[1:nt-1]+y_vec[2:nt]
true_c <- sum(width*height/2, na.rm = T)

save(true_auc_sort, true_c, true_marg_ft, true_marg_ft, file = here("outputData/true_values.RData"))

