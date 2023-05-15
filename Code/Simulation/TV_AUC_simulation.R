
# This script implements the simulation study 
# investigating the effect of overfit

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
set.seed(726)


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

##### simulation ####

# set up: 
# three true covariates
# case 1: fit true model
# case 2: model with 20 additional noise
# case 3: model with 100 additional noise
# case 4: true model, but out-of-sample covariates are generated from different distribution
#         (little overlap)


iter <- 1
skip <- 0 # This is to re-generate data if a dataset has fitting issues
          # because small number of unique events 

pb <- txtProgressBar(min=0, max=M,style=3)
while(iter <= M){
  
  # generate data
  X  <- matrix(rnorm(N*3), ncol=p, nrow=N) # true covariates
  Z <- matrix(rnorm(N*100), ncol = 100, nrow = N) # noise
  data <- gen_St(eta=X %*% Beta, lambda=2, p=2, 
                 gen_Ct = function(N){
                   sample(c(0.5, 1), size = N, replace = T)})
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
  fit_2_lasso_cv <- cv.glmnet(x = cbind(data_train$X, data_train$Z[, 1:20]), 
                        y = Surv(time = data_train$time, event = data_train$event, type = "right"), 
                        family = "cox")
  fit_2_lasso <- glmnet(x = cbind(data_train$X, data_train$Z[, 1:20]), 
                        y= Surv(time = data_train$time, event = data_train$event, type = "right"), 
                        family = "cox", lambda =  fit_2_lasso_cv$lambda.min)
    
  fit_3 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z, data=data_train),
                    warning = function(x){NULL})
  fit_3_lasso_cv <- cv.glmnet(x = cbind(data_train$X, data_train$Z), 
                              y= Surv(time = data_train$time, event = data_train$event, type = "right"), 
                              family = "cox")
  fit_3_lasso <- glmnet(x = cbind(data_train$X, data_train$Z), 
                        y= Surv(time = data_train$time, event = data_train$event, type = "right"), 
                        family = "cox", lambda =  fit_2_lasso_cv$lambda.min)
  
  # if all of the models covereged
  if(!is_null(fit_1) & !is.null(fit_2) & !is.null(fit_3)){
  
    ## estimated in-sample marginal survival
    KM_fit_train <- survfit(Surv(time,event)~1, timefix=FALSE, data=data_train)
    KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]
  
    ## unique time and biomarker values of training sample
    ut_train <- unique(data_train$time[data_train$event==1])
    nt_train <- length(ut_train)
  
    ## in-sample AUC 
    auc_est1 <- tv_auc(eta=fit_1$linear.predictors, data=data_train, t=ut_train, nt=nt_train)
    auc_est2 <- tv_auc(eta=fit_2$linear.predictors, data=data_train, t=ut_train, nt=nt_train)
    auc_est2_lasso <- tv_auc(eta=cbind(data_train$X, data_train$Z[, 1:20]) %*% coef(fit_2_lasso)[,1],
                                data=data_train, t=ut_train, nt=nt_train)
    auc_est3 <- tv_auc(eta=fit_3$linear.predictors, data=data_train, t=ut_train, nt=nt_train)
    auc_est3_lasso <- tv_auc(eta=cbind(data_train$X, data_train$Z) %*% coef(fit_3_lasso)[,1], 
                             data=data_train, t=ut_train, nt=nt_train )
    
    ## in-sample concordance
    concord1 <- concord(data_train, KM_est_train, auc_est1, fit_1$coefficients, data_train$X)
    concord2 <- concord(data_train, KM_est_train, auc_est2, fit_2$coefficients,
                        cbind(data_train$X,data_train$Z[,1:20]))
    concord2_lasso <- concord(data_train, KM_est_train, auc_est2_lasso, coef(fit_2_lasso), 
                              cbind(data_train$X,data_train$Z[,1:20]))
    concord3 <- concord(data_train, KM_est_train, auc_est3, fit_3$coefficients, cbind(data_train$X,data_train$Z))
    concord3_lasso <- concord(data_train, KM_est_train, auc_est3_lasso, coef(fit_3_lasso), 
                              cbind(data_train$X,data_train$Z))
    
    auc_est_train <- rbind(auc_est1, auc_est2, auc_est2_lasso, auc_est3, auc_est3_lasso) %>%
      data.frame() %>%
      mutate(model = rep(c("No noise", "20 noise", "20 noise lasso", "100 noise", "100 noise lasso"), 
                         each = nt_train))
    concord_est_train <- rbind(concord1, concord2, concord2_lasso, concord3, concord3_lasso) %>%
      data.frame() %>%
      mutate(model = c("No noise", "20 noise", "20 noise lasso", "100 noise", "100 noise lasso")) 
  
    # out-of-sample 
    ## estimated risk scores
    ut_test <- unique(data_test$time[data_test$event==1])
    nt_test <- length(ut_test)
  
    ## marginal survival estimates
    KM_fit_test <- survfit(Surv(time,event)~1, timefix=FALSE,data=data_test)
    KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]
  
    ## out-of-sample AUC 
    test_auc_est1 <- tv_auc(data_test$X %*% coef(fit_1), data_test, ut_test, nt_test)
    test_auc_est2 <- tv_auc(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2), data_test, ut_test, nt_test)
    test_auc_est2_lasso <- tv_auc(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2_lasso),
                                  data_test, ut_test, nt_test)
    test_auc_est3 <- tv_auc(cbind(data_test$X,data_test$Z) %*% coef(fit_3), data_test, ut_test, nt_test)
    test_auc_est3_lasso <- tv_auc(cbind(data_test$X,data_test$Z) %*% coef(fit_3_lasso),
                                  data_test, ut_test, nt_test)
  
    ## out-of-sample concordance estimates
    test_concord1 <- concord(data_test, KM_est_test, test_auc_est1, fit_1$coefficients, data_test$X)
    test_concord2 <- concord(data_test, KM_est_test, test_auc_est2, fit_2$coefficients,
                             cbind(data_test$X,data_test$Z[,1:20]))
    test_concord2_lasso <- concord(data_test, KM_est_test, test_auc_est2_lasso, coef(fit_2_lasso), 
                                   cbind(data_test$X,data_test$Z[,1:20]))
    test_concord3 <- concord(data_test, KM_est_test, test_auc_est3,  fit_3$coefficients, cbind(data_test$X,data_test$Z))
    test_concord3_lasso <- concord(data_test, KM_est_test, test_auc_est3_lasso, coef(fit_3_lasso), 
                                   cbind(data_test$X,data_test$Z))
    # clean 
    auc_est_test <- rbind(test_auc_est1, 
                          test_auc_est2,
                          test_auc_est2_lasso, 
                          test_auc_est3, 
                          test_auc_est3_lasso) %>%
      data.frame() %>%
      mutate(model = rep(c("No noise", "20 noise", "20 noise lasso", "100 noise", "100 noise lasso"), 
                       each = nt_test))
    concord_est_test <- rbind(test_concord1, test_concord2, test_concord2_lasso, test_concord3, test_concord3_lasso) %>%
      data.frame() %>%
      mutate(model = c("No noise", "20 noise", "20 noise lasso", "100 noise", "100 noise lasso"))
    
    # put in and out-of-sample together
    tv_auc_df <- bind_rows(
      auc_est_train %>% mutate(sample = "In-sample"), 
      auc_est_test %>% mutate(sample = "Out-of-sample")
    )
    
    c_df <- bind_rows(
      concord_est_train %>% mutate(sample = "In-sample"),
      concord_est_test %>% mutate(sample = "Out-of-sample")
    )
    
    
    # save to final results
    auc_lst[[iter]] <- tv_auc_df
    c_lst[[iter]] <- c_df
    
    # move to next iter
    iter <- iter + 1
  
  }
  # skip non-converged data
  else{skip <- skip+1}

  setTxtProgressBar(pb, value=iter)
}


#### one-iter check #####

save(data_train, data_test, file = here("Data/OneIterData.RData"))


auc_df <- bind_rows(auc_est_train, auc_est_test, .id = "sample")

table(auc_df$model)
head(auc_df)

auc_df_long <- auc_df %>% 
  mutate(# true = true_auc_sort,
    sample = factor(sample, levels = 1:2, labels = c("In-sample", "Out-of-sample"))) %>% 
  pivot_longer(3:5) 


auc_df_long$model <- factor(auc_df_long$model, levels =  c("No noise", "20 noise", "20 noise, penalized", 
                                                           "100 noise", "100 noise, penalized",
                                                           "Different covariate"))
table(auc_df_long$model)
auc_df_long$name <- factor(auc_df_long$name, levels = c("HZ", "SNP", "NP"))
table(auc_df_long$name)


auc_df_long %>% 
  filter(model!="20 noise, penalized" & model!="100 noise, penalized" & model!="Different covariate") %>%
  ggplot(aes(x=time, y=value, col=model, linetype = sample))+
  geom_line()+
  # geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
  #             method = "gam")+
  # geom_line(aes(x = time, y = true), na.rm = T, col = "red", show.legend = F)+
  # labs(x="time", y = "AUC")+
  # theme(text = element_text(size = 15), axis.text = element_text(size = 10))+
  facet_wrap(~name)+
  scale_colour_manual(values=cbPalette)


auc_df_long2 %>% 
  filter(iter == 1) %>% 
  ggplot(aes(x=time, y=value, col=model))+
  geom_point()+
  # geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
  #             method = "gam")+
  # geom_line(aes(x = time, y = true), na.rm = T, col = "red", show.legend = F)+
  # labs(x="time", y = "AUC")+
  # theme(text = element_text(size = 15), axis.text = element_text(size = 10))+
  facet_wrap(~name)+
  scale_colour_manual(values=cbPalette)


table(auc_lst_test_diffx[[1]]$model)

##### All iterations #####

auc_df <- bind_rows(auc_lst, .id="iter")
c_df <- bind_rows(c_lst, .id="iter")

#### true AUC and concordance #####
load(here("Data/true_values.RData"))

# use interpolation to estimate AUC on simulated time series
true_auc_df <- approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, 
                       xout = unique(auc_df$time))
true_auc_df <- data.frame(true_auc_df) %>% 
  rename("time" = x, "true_auc" = y)
head(true_auc_df)

#### Figures ####

# trend of TV-AUC estimates (smoothed)
auc_df %>%
  left_join(true_auc_df, by = "time") %>% 
  pivot_longer(3:5, names_to = "estimator", values_to = "AUC") %>%
  ggplot(aes(x=time, y=AUC, col=model, linetype = sample))+
  geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
              method = "gam")+
  geom_line(aes(x = time, y = true_auc), na.rm = T, col = "red", show.legend = F)+
  labs(x="time", y = "AUC")+
  theme(text = element_text(size = 15), axis.text = element_text(size = 10))+
  facet_wrap(~estimator)+
  scale_colour_manual(values=cbPalette)
ggsave(filename = "tvauc_N250.png", path = "Images/N250", width=15, height=4, bg = "white")
#ggsave(filename = "tvauc_N500.png", path = "Images/N500", width=15, height=4, bg = "white")

# spread of TV-AUC estimates
brk <- seq(0, 1, 0.2) # bin time into five different bins
ggarrange(
auc_df_long %>%
  # dplyr::select(-true) %>%
  filter(sample == "In-sample") %>%
  mutate(time_bin = cut(time, breaks = brk, include.lowest = T)) %>%
  ggplot(aes(x = factor(time_bin), y = value, fill = name))+
  geom_boxplot(outlier.size = 0.5)+
  facet_grid(cols = vars(model))+
  labs(x = "Time", y = "AUC", title = "In-sample")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.1, hjust = 0.1), 
        text=element_text(size = 15),
        axis.text = element_text(size = 10))+
  scale_fill_manual(values=cbPalette)+
  labs(y="AUC"), 

auc_df_long %>%
  # dplyr::select(-true) %>%
  filter(sample == "Out-of-sample") %>%
  mutate(time_bin = cut(time, breaks = brk, include.lowest = T)) %>%
  ggplot(aes(x = factor(time_bin), y = value, fill = name))+
  geom_boxplot(outlier.size = 0.5)+
  facet_grid(cols = vars(model))+
  labs(x = "Time", y = "AUC", title = "Out-of-sample")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.1, hjust = 0.1), 
        text = element_text(size=15),
        axis.text = element_text(size = 10))+
  scale_fill_manual(values=cbPalette)+
  labs(y="AUC"), nrow=1, common.legend = T)
ggsave(filename = "tvauc_box_N250.png", path = "Images/N250", width=15, height = 4, bg = "white")
#ggsave(filename = "tvauc_box_N500.png", path = "Images/N500", width=15, height = 4, bg = "white")

# concordance
## ovefit
c_df %>% pivot_longer(2:6, names_to = "estimator", values_to = "c") %>%
  ggplot()+
  geom_boxplot(aes(x=sample, y=c))+
  facet_grid(rows=vars(estimator), cols=vars(model))
ggsave(filename = "concordance_N250.png", path = "Images/N250", width=15, height = 4, bg = "white")
#ggsave(filename = "concordance_N500.png", path = "Images/N500", width=15, height = 4, bg = "white")


# save(auc_df_train, auc_df_test, c_df_train, c_df_test, file = here("outputData/estimated_values.RData"))
save(auc_lst, c_lst, file = here("Data/SimNoiseModel.RData"))
# load(here("Data/results_by_iter.RData"))


