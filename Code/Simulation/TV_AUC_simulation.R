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
set.seed(71)


# helper functions
# list.files("Code/Simulation")
source(here("Code/Simulation/helpers.R"))
source("Code/Simulation/helpers_estimator.R")

#### Simulation set up ####

M <- 1000 # number of simulated data sets
N <- 500 # number of subjects
Beta <- c(1,-1,0.25)

# training and test sample size 
N_obs <- N*0.5

# number of covariates 
p <- 3

# result container
auc_lst_train <- list()
auc_lst_test <- list()
c_lst_train <- list()
c_lst_test <- list()
data_lst <- list()

##### simulation ####
iter <- 1
skip <- 0

pb <- txtProgressBar(min=0, max=M,style=3)
while(iter <= M){
  
  # generate data
  X  <- matrix(rnorm(N*3), ncol=p, nrow=N)
  Z <- matrix(rnorm(N*100), ncol = 100, nrow = N)
  data <- gen_St(eta=X %*% Beta,lambda=2, p=2, 
                 gen_Ct = function(N){
                   sample(c(0.5, 1), size = N, replace = T)})
  data <- rename(data, true_eta = eta)
  data$X <- I(X)
  data$Z <- I(Z)
  data_lst[[iter]] <- data
  ## separate out the training and test datasets
  data_train <- data[1:N_obs,]
  data_test  <- data[-c(1:N_obs),]
  
  # fit our 3 models with different number of noise predictors (0, 20, 100)
  fit_1 <- tryCatch(expr = coxph(Surv(time, event) ~ X, data=data_train), 
                    warning = function(x){NULL})
  fit_2 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z[,1:20], data=data_train), 
                    warning = function(x){NULL})
  fit_3 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z, data=data_train),
                    warning = function(x){NULL})
  
  # if all of the models covereged
  if(!is_null(fit_1) & !is.null(fit_2) & !is.null(fit_3)){
  
    ## estimated in-sample marginal survival
    KM_fit_train <- survfit(Surv(time,event)~1, timefix=FALSE, data=data_train)
    KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]
  
    ## unique time and biomarker values of training sample
    t_uni_train <- unique(data_train$time[data_train$event==1])
    nt_uni_train <- length(t_uni_train)
  
    ## in-sample AUC and concordance estimates
    auc_est1 <- train_auc(fit_1)
    auc_est2 <- train_auc(fit_2)
    auc_est3 <- train_auc(fit_3)
    concord1 <- train_concord(auc_mat = auc_est1, fit = fit_1)
    concord2 <- train_concord(auc_mat = auc_est2, fit = fit_2)
    concord3 <- train_concord(auc_mat = auc_est3, fit = fit_3)
    auc_est_train <- rbind(auc_est1, auc_est2, auc_est3) %>%
      data.frame() %>%
      mutate(model = rep(c("No noise", "20 noise", "100 noise"), 
                         each = nt_uni_train))
    concord_est_train <- rbind(concord1, concord2, concord3) %>%
      data.frame() %>%
      mutate(model = c("No noise", "20 noise", "100 noise")) 
  
    # out-of-sample 
    ## estimated risk scores
    test_eta1 <- as.vector(data_test$X %*% coef(fit_1))
    test_eta2 <- as.vector(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2))
    test_eta3 <- as.vector(cbind(data_test$X,data_test$Z) %*% coef(fit_3))
    t_uni_test <- unique(data_test$time[data_test$event==1])
    nt_uni_test <- length(t_uni_test)
  
    ## marginal survival estimates
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
  
    # save to final results
    auc_lst_train[[iter]] <- auc_est_train
    auc_lst_test[[iter]] <- auc_est_test
    c_lst_train[[iter]] <- concord_est_train
    c_lst_test[[iter]] <- concord_est_test
    
    # move to next iter
    iter <- iter + 1
  
  }
  # skip non-converged data
  else{skip <- skip+1}

  setTxtProgressBar(pb, value=iter)
}

auc_lst_train[[1]] %>% head()
auc_lst_test[[1]] %>% head()
c_lst_train[[1]]
c_lst_test[[1]]

auc_df_train <- bind_rows(auc_lst_train, .id = "iter")
auc_df_test <- bind_rows(auc_lst_test, .id = "iter")
auc_df <- bind_rows(auc_df_train, auc_df_test, .id = "sample")
unique(auc_df_train$iter)


## concordance
c_df_train <- bind_rows(c_lst_train, .id = "iter")
c_df_test <- bind_rows(c_lst_test, .id = "iter")
c_df <- bind_rows(c_df_train, c_df_test, .id = "sample")

#### true AUC #####
load(here("outputData/true_values.RData"))
true_auc_sort <- approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, 
                       xout = auc_df$time)$y

head(true_auc_sort)
plot(auc_df$time, true_auc_sort)

# interpolated_auc_train<-approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, 
#                                xout = auc_df_train$time)$y
# interpolated_auc_test<-approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, 
#                               xout = auc_df_test$time)$y

#### brief check of results #####
auc_df %>% 
  mutate(true = true_auc_sort) %>% 
  pivot_longer(4:6) %>% 
  ggplot(aes(x=time, y=value, col=model, linetype = sample))+
  geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
              method = "gam")+
  facet_wrap(~name)+
  geom_line(aes(x=time, y = true))
  
  
c_df <- bind_rows(c_df_train, c_df_test, .id = "sample")
c_df %>%
  pivot_longer(3:7) %>%
  ggplot(aes(x=name, y=value))+
  geom_boxplot()+
  facet_wrap(~sample)

save(auc_df, c_df, file = here("outputData/estimates_N_250.RData"))

# save(auc_df_train, auc_df_test, c_df_train, c_df_test, file = here("outputData/estimated_values.RData"))
# save(auc_lst_train, auc_lst_test, file = here("outputData/results_by_iter.RData"))

