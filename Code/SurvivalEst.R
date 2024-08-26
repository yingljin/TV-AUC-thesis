# This script explores the influnce of basis function
# on the smoothed survival function

#### Set up ####

rm(list = ls())
library("survival")
library("ggplot2")
library(ggpubr)
library("cubature")
library(tidyverse)
library(dplyr)
library(here)
library(risksetROC)
library(mgcv)
library(scam)
library("clinfun")
library(glmnet)
library(kableExtra)
theme_set(theme_minimal())
cbPalette <- c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set the seed
set.seed(512)

load(here("Data/SimNoiseTrainData.RData"))

df <- train_list[[1]] 

#### Survival curves ####

# KP
kp_surv <- survfit(Surv(time,event)~1, timefix=FALSE, data=df)
table(df$event)
plot(kp_surv)
kp_surv$surv[kp_surv$n.event>0]

df_surv <- data.frame(time = kp_surv$time[kp_surv$n.event>0],
                      kp = kp_surv$surv[kp_surv$n.event>0])

df_surv2 <- expand_grid(df_surv, 
                       k = c(10, 20, 30),
                       basis = c("MPD", "MDCX"))

# SCAM smoothing
# shape constrained smoothing terms (SCOP-splines)

k_vec <- c(10, 20, 30)
df_surv2$SmSurv <- NA

for(i in seq_along(k_vec)){
  # MPD
  mpd_surv <- scam(kp ~ s(time, bs="mpd",k=k_vec[i]),
                    data=df_surv)
  df_surv2$SmSurv[df_surv2$basis=="MPD"&df_surv2$k==k_vec[i]] <- predict(mpd_surv)
  
  # MPCX
  mpcx_surv <- scam(kp ~ s(time, bs="mdcx",k=k_vec[i]),
                   data=df_surv)
  df_surv2$SmSurv[df_surv2$basis=="MDCX"&df_surv2$k==k_vec[i]] <- predict(mpcx_surv)
  
  
}



# Figure

df_surv2 %>% 
  mutate(k = paste0("k = ", k)) %>%
  ggplot()+
  geom_step(aes(x=time, y=kp, col = "Kaplan-Meier"))+
  geom_line(aes(x=time, y=SmSurv, col = "Smooth"))+
  labs(y = "", x="Time", col = "")+
  scale_color_manual(values = cbPalette)+
  facet_grid(rows = vars(basis), cols = vars(k))
ggsave(filename = here("Code/SurvSmooth.pdf"), 
       width=12, height=8, bg="white", dpi = 300)


#### How basis affects concordance ####

# Review round 1, review 2's question 2

##### first, calculated I/D AUC #####

# helper functions
source(here("Code/Simulation/helpers.R"))

# redefine concordance function

# set up
M <- 1000 # number of simulated data sets
N <- 500 # number of subjects
Beta <- c(1,-1,0.25) # true coefficients

# training and test sample size 
N_obs <- N*0.5

# number of covariates 
p <- 3



# result container
auc_lst <- list()
# c_lst <- list()
# train_list <- list()
# test_list <- list()

# set up: 
# three true covariates
# case 1: fit true model
# case 2: model with 20 additional noise
#         PS: lasso penalization was also investigated
# case 3: model with 100 additional noise
#         PS: lasso penalization, was also investigated

iter <- 1
skip <- 0 # This is to re-generate data if a dataset has fitting issues
# because small number of unique events 

t1 <- Sys.time()
pb <- txtProgressBar(min=0, max=M,style=3)
while(iter <= M){
  
  # generate data
  X  <- matrix(rnorm(N*3), ncol=p, nrow=N) # true covariates
  Z <- matrix(rnorm(N*100), ncol = 100, nrow = N) # noise
  data <- gen_St(eta=X %*% Beta, lambda=2, p=2, 
                 gen_Ct = function(N){runif(N, 0, 1)})
  # gen_Ct = function(N){sample(c(0.5, 1), size = N, replace = T)}
  data <- rename(data, true_eta = eta)
  data$X <- I(X)
  data$Z <- I(Z)
  # data_lst[[iter]] <- data
  ## separate out the training and test datasets
  data_train <- data[1:N_obs,]
  data_test  <- data[-c(1:N_obs),]
  
  # fit our 3 models with different number of noise predictors (0, 20, 100)
  # and the latter two with lasso penalization
  fit_1 <- tryCatch(expr = coxph(Surv(time, event) ~ X, data=data_train), 
                    warning = function(x){NULL})
  
  fit_2 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z[,1:20], data=data_train), 
                    warning = function(x){NULL})
  # fit_2_lasso_cv <- cv.glmnet(x = cbind(data_train$X, data_train$Z[, 1:20]), 
  #                       y = Surv(time = data_train$time, event = data_train$event, type = "right"), 
  #                       family = "cox")
  # fit_2_lasso <- glmnet(x = cbind(data_train$X, data_train$Z[, 1:20]), 
  #                       y= Surv(time = data_train$time, event = data_train$event, type = "right"), 
  #                       family = "cox", lambda =  fit_2_lasso_cv$lambda.min)
  
  fit_3 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z, data=data_train),
                    warning = function(x){NULL})
  # fit_3_lasso_cv <- cv.glmnet(x = cbind(data_train$X, data_train$Z), 
  #                             y= Surv(time = data_train$time, event = data_train$event, type = "right"), 
  #                             family = "cox")
  # fit_3_lasso <- glmnet(x = cbind(data_train$X, data_train$Z), 
  #                       y= Surv(time = data_train$time, event = data_train$event, type = "right"), 
  #                       family = "cox", lambda =  fit_2_lasso_cv$lambda.min)
  
  # if all of the models covereged
  if(!is_null(fit_1) & !is.null(fit_2) & !is.null(fit_3)){
    
    ## estimated in-sample marginal survival
    KM_fit_train <- survfit(Surv(time,event)~1, timefix=FALSE, data=data_train)
    # KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]
    
    ## unique time and biomarker values of training sample
    ut_train <- unique(data_train$time[data_train$event==1])
    nt_train <- length(ut_train)
    
    ## in-sample AUC 
    auc_est1 <- tv_auc(eta=fit_1$linear.predictors, data=data_train, t=ut_train, nt=nt_train)
    auc_est2 <- tv_auc(eta=fit_2$linear.predictors, data=data_train, t=ut_train, nt=nt_train)
    auc_est3 <- tv_auc(eta=fit_3$linear.predictors, data=data_train, t=ut_train, nt=nt_train)
    
    ## add in survival probability
    KM_est_train <- data.frame(time = KM_fit_train$time[KM_fit_train$n.event>0],
                               surv = KM_fit_train$surv[KM_fit_train$n.event>0]) 
    
    auc_est_train <- rbind(auc_est1, auc_est2, auc_est3) %>%
      data.frame() %>%
      mutate(model = rep(c("No noise", "20 noise", "100 noise"), 
                         each = nt_train)) %>% 
      left_join(KM_est_train, by = "time")
    
    
    # out-of-sample 
    ## estimated risk scores
    ut_test <- unique(data_test$time[data_test$event==1])
    nt_test <- length(ut_test)
    
    ## marginal survival estimates
    KM_fit_test <- survfit(Surv(time,event)~1, timefix=FALSE,data=data_test)
    # KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]
    
    ## out-of-sample AUC 
    test_auc_est1 <- tv_auc(data_test$X %*% coef(fit_1), data_test, ut_test, nt_test)
    test_auc_est2 <- tv_auc(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2), data_test, ut_test, nt_test)
    test_auc_est3 <- tv_auc(cbind(data_test$X,data_test$Z) %*% coef(fit_3), data_test, ut_test, nt_test)
    
    ## add in survival probability
    KM_est_test <- data.frame(time = KM_fit_test$time[KM_fit_test$n.event>0],
                               surv = KM_fit_test$surv[KM_fit_test$n.event>0]) 
    
    # clean 
    auc_est_test <- rbind(test_auc_est1, 
                          test_auc_est2,
                          # test_auc_est2_lasso, 
                          test_auc_est3) %>%
      # test_auc_est3_lasso) %>%
      data.frame() %>%
      mutate(model = rep(c("No noise", "20 noise","100 noise"), 
                         each = nt_test)) %>%
      left_join(KM_est_test, by = "time")
    # concord_est_test <- rbind(test_concord1, test_concord2,test_concord3) %>%
    #   data.frame() %>%
    #   mutate(model = c("No noise", "20 noise", "100 noise"))
    
    # put in and out-of-sample together
    tv_auc_df <- bind_rows(
      auc_est_train %>% mutate(sample = "In-sample"), 
      auc_est_test %>% mutate(sample = "Out-of-sample")
    )
    
    # c_df <- bind_rows(
    #   concord_est_train %>% mutate(sample = "In-sample"),
    #   concord_est_test %>% mutate(sample = "Out-of-sample")
    # )
    # 
    
    # save to final results
    auc_lst[[iter]] <- tv_auc_df
    # c_lst[[iter]] <- c_df
    # train_list[[iter]] <- data_train
    # test_list[[iter]] <- data_test
    
    # move to next iter
    iter <- iter + 1
    
  }
  # skip non-converged data
  else{skip <- skip+1}
  
  setTxtProgressBar(pb, value=iter)
}


close(pb)
t2 <- Sys.time()

##### - second, integrate AUC #####

# define a function that can use different splines
intAUC_smS <- function(AUC, ut, St, k=10, bs_type){
  # smooth survivla function
  fit_S <- scam(St ~ s(ut, bs=bs_type,k=k), data=data.frame(ut=ut, St=St))
  df_pred <- data.frame(ut = rep(ut, each=2) + rep(c(-0.00001,0.00001), length(ut)))
  St_pred <- predict(fit_S, newdata=df_pred, type='response')
  ft     <- -diff(St_pred)[seq(1,2*length(ut),by=2)]                    
  
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

View(auc_lst[[1]])
auc_lst[[1]] %>% arrange(sample, model, time) %>% 
  group_by(sample, model) %>%
  group_modify(~{
    data.frame(C = intAUC_smS(.x$SNP, .x$time, .x$surv, k=10, bs_type = "mdcx"))
  })

auc_sort <-arrange(data.frame(auc_lst[[1]] %>% ), time)
intAUC(auc_sort$SNP, auc_sort$time, KM_est, 
       method = "smS") 

# weighted by survival functions smoothed by different splines

NP_C_list <- list()
SNP_C_list <- list()

pb <- txtProgressBar(min=0, max=M,style=3)

t1 <- Sys.time()
## non-parametric
for(m in seq_along(auc_lst)){
  
  NP_C_df <- tryCatch(
    expr = {auc_lst[[m]] %>% arrange(sample, model, time) %>% 
              group_by(sample, model) %>%
              group_modify(~{
                data.frame(
                  # MPD
                  C_k10_mpd_np = intAUC_smS(.x$NP, .x$time, .x$surv, k=10, bs_type = "mpd"),
                  C_k20_mpd_np = intAUC_smS(.x$NP, .x$time, .x$surv, k=20, bs_type = "mpd"),
                  C_k30_mpd_np = intAUC_smS(.x$NP, .x$time, .x$surv, k=30, bs_type = "mpd"),
                  # MDCX
                  C_k10_mdcx_np = intAUC_smS(.x$NP, .x$time, .x$surv, k=10, bs_type = "mdcx"),
                  C_k20_mdcx_np = intAUC_smS(.x$NP, .x$time, .x$surv, k=20, bs_type = "mdcx"),
                  C_k30_mdcx_np = intAUC_smS(.x$NP, .x$time, .x$surv, k=30, bs_type = "mdcx"))
              })},
    error = function(e){NA}
    
    )
  
  SNP_C_df <- tryCatch(
    expr = {auc_lst[[m]] %>% arrange(sample, model, time) %>% 
              group_by(sample, model) %>%
              group_modify(~{
                data.frame(
                  # MPD
                  C_k10_mpd_np = intAUC_smS(.x$SNP, .x$time, .x$surv, k=10, bs_type = "mpd"),
                  C_k20_mpd_np = intAUC_smS(.x$SNP, .x$time, .x$surv, k=20, bs_type = "mpd"),
                  C_k30_mpd_np = intAUC_smS(.x$SNP, .x$time, .x$surv, k=30, bs_type = "mpd"),
                  # MDCX
                  C_k10_mdcx_np = intAUC_smS(.x$SNP, .x$time, .x$surv, k=10, bs_type = "mdcx"),
                  C_k20_mdcx_np = intAUC_smS(.x$SNP, .x$time, .x$surv, k=20, bs_type = "mdcx"),
                  C_k30_mdcx_np = intAUC_smS(.x$SNP, .x$time, .x$surv, k=30, bs_type = "mdcx"))
              })},
    error = function(e){NA})
  
  NP_C_list[[m]] <- NP_C_df
  SNP_C_list[[m]] <- SNP_C_df
  
  setTxtProgressBar(pb, value=m)
  
}

close(pb)
t2 <- Sys.time()

##### results ####
# remove the error ones
NP_C_list <- NP_C_list[!is.na(NP_C_list)]
SNP_C_list <- SNP_C_list[!is.na(SNP_C_list)]

bind_rows(
  bind_rows(NP_C_list, .id = "iter") %>% mutate(est = "Non-parametric"),
  bind_rows(SNP_C_list, .id = "iter") %>% mutate(est = "Smoothed non-parametric")
) %>%
  pivot_longer(starts_with("C_k")) %>%
  ggplot()+
  geom_boxplot(aes(x=model, y = value, fill = name))+
  facet_grid(rows = vars(est), cols = vars(sample))+
  labs(y = "Concordance", x = "", fill = "Basis")+
  scale_fill_manual(values = cbPalette,
                    labels = c("K=10, MDCX", "K=10, MPD", 
                               "K=20, MDCX", "K=20, MPD", 
                               "K=30, MDCX", "K=30, MPD"))

ggsave(filename = here("Code/ConcordSurvSmooth.pdf"), 
       width=10, height=8, bg="white", dpi = 300)
