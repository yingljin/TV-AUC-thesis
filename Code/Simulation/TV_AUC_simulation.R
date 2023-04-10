
# This script implements the simulation study that 
# corresponding to the simulation section 
# in the manuscript

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
auc_lst_train <- list()
auc_lst_test <- list()
c_lst_train <- list()
c_lst_test <- list()
data_lst <- list()

##### simulation ####
iter <- 1
skip <- 0 # This is to re-generate data if a dataset has fitting issues
          # because small number of unique events 

pb <- txtProgressBar(min=0, max=M,style=3)
while(iter <= M){
  
  # generate data
  X  <- matrix(rnorm(N*3), ncol=p, nrow=N)
  Z <- matrix(rnorm(N*100), ncol = 100, nrow = N)
  data <- gen_St(eta=X %*% Beta, lambda=2, p=2, 
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
  # and the latter two with penalization
  fit_1 <- tryCatch(expr = coxph(Surv(time, event) ~ X, data=data_train), 
                    warning = function(x){NULL})
  
  fit_2 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z[,1:20], data=data_train), 
                    warning = function(x){NULL})
  
  # cross validation paritial likelihood -based criterion
  fit_2_lasso_cv <- cv.glmnet(x = cbind(data_train$X, data_train$Z[, 1:20]), 
                        y= Surv(time = data_train$time, event = data_train$event, type = "right"), 
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
    t_uni_train <- unique(data_train$time[data_train$event==1])
    nt_uni_train <- length(t_uni_train)
  
    ## in-sample AUC and concordance estimates
    auc_est1 <- train_auc(eta=fit_1$linear.predictors)
    auc_est2 <- train_auc(eta=fit_2$linear.predictors)
    auc_est2_lasso <- train_auc(eta=cbind(data_train$X, data_train$Z[, 1:20]) %*% coef(fit_2_lasso)[,1])
    auc_est3 <- train_auc(eta=fit_3$linear.predictors)
    auc_est3_lasso <- train_auc(eta=cbind(data_train$X, data_train$Z) %*% coef(fit_3_lasso)[,1] )
    
    
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
    test_eta2_lasso <- as.vector(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2_lasso))
    test_eta3 <- as.vector(cbind(data_test$X,data_test$Z) %*% coef(fit_3))
    test_eta3_lasso <- as.vector(cbind(data_test$X,data_test$Z) %*% coef(fit_3_lasso))
    t_uni_test <- unique(data_test$time[data_test$event==1])
    nt_uni_test <- length(t_uni_test)
  
    ## marginal survival estimates
    KM_fit_test <- survfit(Surv(time,event)~1, timefix=FALSE,data=data_test)
    KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]
  
    ## in-sample AUC and concordance estimates
    test_auc_est1 <- test_auc(test_eta1)
    test_auc_est2 <- test_auc(test_eta2)
    test_auc_est2_lasso <- test_auc(test_eta2_lasso)
    test_auc_est3 <- test_auc(test_eta3)
    test_auc_est3_lasso <- test_auc(test_eta3_lasso)
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

##### Cleaning and formatting results #####

## AUC
auc_df_train <- bind_rows(auc_lst_train, .id = "iter")
auc_df_test <- bind_rows(auc_lst_test, .id = "iter")
auc_df <- bind_rows(auc_df_train, auc_df_test, .id = "sample")


## concordance
c_df_train <- bind_rows(c_lst_train, .id = "iter")
c_df_test <- bind_rows(c_lst_test, .id = "iter")
c_df <- bind_rows(c_df_train, c_df_test, .id = "sample")

#### true AUC and concordance #####
load(here("outputData/true_values.RData"))

# use interpolation to estimate AUC on simulated time series
true_auc_sort <- approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, 
                       xout = auc_df$time)$y



#### Produce figures #####

# trend of TV-AUC estimates (smoothed)
auc_df_long <- auc_df %>% 
  mutate(true = true_auc_sort,
         sample = factor(sample, levels = 1:2, labels = c("In-sample", "Out-of-sample"))) %>% 
  pivot_longer(4:6) 

auc_df_long$model <- factor(auc_df_long$model, levels = c("No noise", "20 noise", "100 noise"))
table(auc_df_long$model)
auc_df_long$name <- factor(auc_df_long$name, levels = c("HZ", "SNP", "NP"))
table(auc_df_long$name)


auc_df_long %>% 
  ggplot(aes(x=time, y=value, col=model, linetype = sample))+
  geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
              method = "gam")+
  geom_line(aes(x = time, y = true), na.rm = T, col = "red", show.legend = F)+
  labs(x="time", y = "AUC")+
  theme(text = element_text(size = 15), axis.text = element_text(size = 10))+
  facet_wrap(~name)+
  scale_colour_manual(values=cbPalette)

ggsave(filename = "tvauc_N250.png", path = "Images/N250", width=15, height=4, bg = "white")
#ggsave(filename = "tvauc_N500.png", path = "Images/N500", width=15, height=4, bg = "white")

# spread of TV-AUC estimates
brk <- seq(0, 1, 0.2) # bin time into five different bins
ggarrange(
auc_df_long %>%
  dplyr::select(-true) %>%
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
  dplyr::select(-true) %>%
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
c_df <- bind_rows(c_df_train, c_df_test, .id = "sample")
c_df_long <- c_df %>%
  mutate(sample = factor(sample, levels = 1:2, labels = c("In-sample", "Out-of-sample"))) %>% 
  pivot_longer(3:7) %>%
  mutate(Type = ifelse(name=="HZ"|name=="GH", "Semi-parametrc", "Non-parametric")) %>%
  mutate(model = factor(model, levels = c("No noise",  "20 noise", "100 noise")),
         name = factor(name, levels = c("HZ","GH", "Harrell.s","NP","SNP"),
                       labels = c("HZ","GH", "Harrell","NP","SNP")))

ggarrange(
c_df_long %>% 
  filter(sample == "In-sample") %>%
  ggplot(aes(x = name, y = value))+
  geom_boxplot(aes(fill = Type))+
  facet_grid(cols=vars(model))+
  geom_hline(yintercept = true_c, col = "red")+
  theme(text = element_text(size = 12),
        axis.text = element_text(size=8))+
  scale_fill_manual(values=cbPalette)+
  labs(y = "Concordance", x = "Estimator", title = "In-sample"),

c_df_long %>% 
  filter(sample == "Out-of-sample") %>%
  ggplot(aes(x = name, y = value))+
  geom_boxplot(aes(fill = Type))+
  facet_grid(cols=vars(model))+
  geom_hline(yintercept = true_c, col = "red")+
  theme(text = element_text(size = 12),
        axis.text = element_text(size=8))+
  scale_fill_manual(values=cbPalette)+
  labs(y = "Concordance", x = "Estimator", title = "Out-of-sample"), nrow = 1, common.legend = T)
ggsave(filename = "concordance_N250.png", path = "Images/N250", width=15, height = 4, bg = "white")
#ggsave(filename = "concordance_N500.png", path = "Images/N500", width=15, height = 4, bg = "white")



# save(auc_df_train, auc_df_test, c_df_train, c_df_test, file = here("outputData/estimated_values.RData"))
# save(auc_lst_train, auc_lst_test, file = here("outputData/results_by_iter.RData"))

