# This script implements one iteration of the simulations study
# with a brief exploration of the simulated data


##### set up #####
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
set.seed(726)

# helper functions
# list.files("Code/Simulation")
source(here("Code/Simulation/helpers.R"))
source("Code/Simulation/helpers_estimator.R")

#### Simulation set up ####

N <- 500 # number of subjects
Beta <- c(1,-1,0.25) # true coefficient

# training and test sample size 
N_obs <- N*0.5

# number of covariates 
p <- 3

# result container
auc_lst_train <- list()
auc_lst_test <- list()
c_lst_train <- list()
c_lst_test <- list()


##### generate data  ####
X <- matrix(rnorm(N*3), ncol=p, nrow=N) # covariates
Z <- matrix(rnorm(N*100), ncol = 100, nrow = N) # noise
data <- gen_St(eta=X %*% Beta, lambda=2, p=2, 
                gen_Ct = function(N){
                  sample(c(0.5, 1), size = N, replace = T)}) # survival and censor
data <- rename(data, true_eta = eta)
data$X <- I(X)
data$Z <- I(Z)

# Kaplan-Meier curve
plot(survfit(Surv(time, event)~1, data = data))

# separate out the training and test datasets
data_train <- data[1:N_obs,]
data_test  <- data[-c(1:N_obs),]
  
#### fit cox models on training data ####
# true model
fit_1 <- tryCatch(expr = coxph(Surv(time, event) ~ X, data=data_train), 
                    warning = function(x){NULL})
plot(fit_1)
# moderately overfitted
fit_2 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z[,1:20], data=data_train), 
                    warning = function(x){NULL})
# severely overfitted  
fit_3 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z, data=data_train),
                    warning = function(x){NULL})
  
#### marginal survival on training data ####
KM_fit_train <- survfit(Surv(time,event)~1, timefix=FALSE, data=data_train)
KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]
    
# unique time and biomarker values of training sample
t_uni_train <- unique(data_train$time[data_train$event==1])
nt_uni_train <- length(t_uni_train)

#### in-sample AUC and concordance estimates  ####
# AUC and concordance estimates
auc_est1 <- train_auc(fit_1)
auc_est2 <- train_auc(fit_2)
auc_est3 <- train_auc(fit_3)
concord1 <- train_concord(auc_mat = auc_est1, fit = fit_1)
concord2 <- train_concord(auc_mat = auc_est2, fit = fit_2)
concord3 <- train_concord(auc_mat = auc_est3, fit = fit_3)

# format results
auc_est_train <- rbind(auc_est1, auc_est2, auc_est3) %>%
  data.frame() %>%
  mutate(model = rep(c("No noise", "20 noise", "100 noise"), 
                     each = nt_uni_train))
concord_est_train <- rbind(concord1, concord2, concord3) %>%
  data.frame() %>%
      mutate(model = c("No noise", "20 noise", "100 noise")) 

# brief overview
ggplot(auc_est_train)+
  geom_line(aes(x = time, y = HZ, col = model))
ggplot(auc_est_train)+
  geom_line(aes(x = time, y = NP, col = model))
ggplot(auc_est_train)+
  geom_line(aes(x = time, y = SNP, col = model))

#### out-of-sample estimates ####
# estimated risk scores
## based on model trained in-sample
test_eta1 <- as.vector(data_test$X %*% coef(fit_1))
test_eta2 <- as.vector(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2))
test_eta3 <- as.vector(cbind(data_test$X,data_test$Z) %*% coef(fit_3))

# unique event time in the test set
t_uni_test <- unique(data_test$time[data_test$event==1])
nt_uni_test <- length(t_uni_test)

# marginal survival estimates on the test set
KM_fit_test <- survfit(Surv(time,event)~1, timefix=FALSE,data=data_test)
KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]

## in-sample AUC and concordance estimates
test_auc_est1 <- test_auc(test_eta1)
test_auc_est2 <- test_auc(test_eta2)
test_auc_est3 <- test_auc(test_eta3)
test_concord1 <- test_concord(auc_mat = test_auc_est1, fit = fit_1, 
                              eta = test_eta1, 
                              X_mat = data_test$X)
test_concord2 <- test_concord(auc_mat = test_auc_est2, fit = fit_2, 
                              eta = test_eta2,
                              X_mat = cbind(data_test$X,data_test$Z[,1:20]))
test_concord3 <- test_concord(auc_mat = test_auc_est3, fit = fit_3, 
                              eta = test_eta3,
                              X_mat = cbind(data_test$X,data_test$Z))
auc_est_test <- rbind(test_auc_est1, test_auc_est2, test_auc_est3) %>%
  data.frame() %>%
  mutate(model = rep(c("No noise", "20 noise", "100 noise"), 
                     each = nt_uni_test))
concord_est_test <- rbind(test_concord1, test_concord2, test_concord3) %>%
  data.frame() %>%
  mutate(model = c("No noise", "20 noise", "100 noise")) 

# brief overview
ggplot(auc_est_test)+
  geom_line(aes(x = time, y = HZ, col = model))
ggplot(auc_est_test)+
  geom_line(aes(x = time, y = NP, col = model))
ggplot(auc_est_test)+
  geom_line(aes(x = time, y = SNP, col = model))


