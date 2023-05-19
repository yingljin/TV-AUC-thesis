# This script implements the simulation study 
# investigating the effect of contamination
# by introducing outliying observations in the testing dataset

##### set up #####
rm(list = ls())
library("survival")
library("ggplot2")
library(ggpubr)
library("cubature")
library("tidyverse")
library(here)
library(risksetROC)
library(mgcv)
library(scam)
library("clinfun")
library(glmnet)
theme_set(theme_minimal())
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# set the seed
set.seed(512)


# helper functions
# list.files("Code/Simulation")
source(here("Code/Simulation/helpers.R"))
source("Code/Simulation/helpers_estimator.R")

#### Simulation set up ####

M <- 1000 # number of simulated data sets
N <- 500 # number of subjects
#N <- 1000 # number of subjects
Beta <- c(1,-1,0.25)

# training and test sample size 
N_obs <- N*0.5

# number of covariates 
p <- 3

# result container
auc_lst <- list()
c_lst <- list()
data_lst <- list()
##### simulation ####

# set up: 
# training data: three true covariates
# testing data: introduce outliers
# case 1: 10% outliers, covariate from N(10, 0)
# case 2: 10% outliers, covariate from N(0, 10)



iter <- 1
# skip <- 0 # This is to re-generate data if a dataset has fitting issues
# because small number of unique events 

pb <- txtProgressBar(min=0, max=M,style=3)
while(iter <= M){
  
  # generate data
  X  <- matrix(rnorm(N*3), ncol=p, nrow=N) # true covariates
  data <- gen_St(eta=X %*% Beta, lambda=2, p=2, 
                 gen_Ct = function(N){
                   sample(c(0.5, 1), size = N, replace = T)})
  data <- rename(data, true_eta = eta)
  data$X <- I(X)
  data_lst[[iter]] <- data
  
  # test data 0: no outlier
  data_train <- data[1:N_obs,]
  data_test  <- data[-c(1:N_obs),]
  
  # test data 1: 10% outlier, 10 times larger risk score
  out_id1 <- sample(1:N_obs, size = 0.1*N_obs)
  data_test1 <- data_test
  data_test1$X[out_id1, ] <- rnorm(0.1*N_obs*3, mean=5, sd=1) # adjust covariate value
  data_test1$true_eta[out_id1] <- data_test1$X[out_id1, ] %*% Beta # risk score
  # boxplot(data_test1$true_eta)
  
  # test data 2: 10% outlier, 100 times larger risk score
  data_test2 <- data_test
  data_test2$X[out_id1, ] <- rnorm(0.1*N_obs*3, mean=0, sd=5) # adjust covariate value
  data_test2$true_eta[out_id1] <- data_test2$X[out_id1, ] %*% Beta # risk score
  # boxplot(data_test2$true_eta)
  
  # test data 3: 20% outlier, 10 times larger risk score
  # out_id2 <- sample(1:N_obs, size = 0.2*N_obs)
  # data_test3 <- data_test
  # data_test3$X[out_id2, ] <- 10*data_test3$X[out_id2, ] # adjust covariate value
  # data_test3$true_eta[out_id2] <- data_test3$X[out_id2, ] %*% Beta # risk score
  #boxplot(data_test3$true_eta)
  
  # fit our cox model on training data
  fit_cox <- coxph(Surv(time, event) ~ X, data=data_train)
  
  # in-sample tv-auc and concordance
  ## in-sample suvival functions
  KM_fit_train <- survfit(Surv(time,event)~1, timefix=FALSE, data=data_train)
  KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]
    
  ## unique time and biomarker values of training sample
  ut_train <- unique(data_train$time[data_train$event==1])
  nt_train <- length(ut_train)
    
  ## in-sample AUC 
  auc_est_train <- tv_auc(eta=fit_cox$linear.predictors, 
                          data=data_train,
                          t=ut_train, nt=nt_train)
  ## in-sample concordance
  c_train <- concord(data_train, KM_est_train, auc_est_train, 
                     fit_cox$coefficients, data_train$X)
  
  # out of sample
  # since only covariates are changed
  # marginal survival functions and unique event time stay unchanged
  # across all cases of test data
  KM_fit_test <- survfit(Surv(time,event)~1, timefix=FALSE,data=data_test)
  KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]
  ut_test <- unique(data_test$time[data_test$event==1])
  nt_test <- length(ut_test)
  
  # out-of-sample TV-AUC and concordance
  ## case 0
  auc_est_test <- tv_auc(data_test$X %*% coef(fit_cox), data_test, 
                         ut_test, nt_test)
  c_test <- concord(data_test, KM_est_test, auc_est_test, 
                    fit_cox$coefficients, data_test$X)
  
  # case 1
  auc_est_test1 <- tv_auc(data_test1$X %*% coef(fit_cox), data_test1, 
                         ut_test, nt_test)
  c_test1 <- concord(data_test1, KM_est_test, auc_est_test1, 
                    fit_cox$coefficients, data_test1$X)
  # case 2
  auc_est_test2 <- tv_auc(data_test2$X %*% coef(fit_cox), data_test2, 
                          ut_test, nt_test)
  c_test2 <- concord(data_test2, KM_est_test, auc_est_test2, 
                     fit_cox$coefficients, data_test2$X)
  
  # case 3
  # auc_est_test3 <- tv_auc(data_test3$X %*% coef(fit_cox), data_test3, 
  #                         ut_test, nt_test)
  # c_test3 <- concord(data_test3, KM_est_test, auc_est_test3, 
  #                    fit_cox$coefficients, data_test3$X)
  
  # clean results
  tv_auc_df <- bind_rows(
    auc_est_train %>% data.frame() %>% mutate(sample = "In-sample"),
    auc_est_test %>% data.frame() %>% mutate(sample = "Out-of-sample", 
                                            outlier = "None"),
    auc_est_test1 %>% data.frame() %>% mutate(sample = "Out-of-sample", 
                                             outlier = "N(5, 1)"),
    auc_est_test2 %>% data.frame() %>% mutate(sample = "Out-of-sample", 
                                             outlier = "N(0, 5)"),
    # auc_est_test3 %>% data.frame() %>% mutate(sample = "Out-of-sample", 
    #                                          outlier = "20%, sd = 10")
  )
  
  c_df <- bind_rows(
    c_train %>% data.frame() %>% mutate(sample = "In-sample"),
    c_test %>%  data.frame() %>% mutate(sample = "Out-of-sample",
                                        outlier = "None"),
    c_test1 %>% data.frame() %>% mutate(sample = "Out-of-sample", 
                                              outlier = "N(5, 1)"),
    c_test2 %>% data.frame() %>% mutate(sample = "Out-of-sample", 
                                              outlier = "N(0, 5)"),
    # c_test3 %>% data.frame() %>% mutate(sample = "Out-of-sample", 
    #                                           outlier = "20%, sd = 10")
  )
    
    # save to final results
    auc_lst[[iter]] <- tv_auc_df
    c_lst[[iter]] <- c_df
    
    
    setTxtProgressBar(pb, value=iter)
    iter <- iter+1
}

##### One iteration #####

save(data_test, data_test3, data_test2, data_test1, 
     file = here("Data/OneIterContamData.RData"))

#### Results ####

save(auc_lst, c_lst, file = here("Data/SimContamData.RData"))
