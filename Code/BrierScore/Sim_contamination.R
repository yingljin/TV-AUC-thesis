# This script implements the simulation study 
# investigating the effect of contamination
# corresponding to Section 4.2 in the manuscript

##### set up #####
library("SurvMetrics")
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
set.seed(513)


# helper functions
# list.files("Code/Simulation")
source(here("Code/Simulation/helpers.R"))

#### Simulation set up ####

M <- 1000 # number of simulated data sets
N <- 500 # number of subjects
Beta <- c(1,-1,0.25) # true coefficients

# training and test sample size 
N_obs <- N*0.5

# number of covariates 
p <- 3

# result container
bs_lst <- list()
# c_lst <- list()
# train_lst <- list()
# test_lst <- list()

##### simulation ####

# set up: 
# training data: three true covariates
# testing data: introduce outliers
# case 1: 10% outliers, covariate from N(5, 0)
# case 2: 10% outliers, covariate from N(0, 5)
# case 3: no outliers

iter <- 1

pb <- txtProgressBar(min=0, max=M, style=3)
t1 <- Sys.time()
while(iter <= M){
  
  # generate data
  X  <- matrix(rnorm(N*3), ncol=p, nrow=N) # true covariates
  data <- gen_St(eta=X %*% Beta, lambda=2, p=2, 
                 gen_Ct = function(N){runif(N, 0, 1)})
  data <- rename(data, true_eta = eta)
  data$X <- I(X)
  # data_lst[[iter]] <- data
  
  # training data
  data_train <- data[1:N_obs,]
  
  # test data 0: no outlier
  data_test  <- data[-c(1:N_obs),]
  
  # test data 1: 10% outlier, mean change
  out_id1 <- sample(1:N_obs, size = 0.1*N_obs)
  data_test1 <- data_test
  data_test1$X[out_id1, ] <- rnorm(0.1*N_obs*3, mean=5, sd=1) # adjust covariate value
  data_test1$true_eta[out_id1] <- data_test1$X[out_id1, ] %*% Beta # risk score
  
  # test data 2: 10% outlier, sd change
  data_test2 <- data_test
  data_test2$X[out_id1, ] <- rnorm(0.1*N_obs*3, mean=0, sd=5) # adjust covariate value
  data_test2$true_eta[out_id1] <- data_test2$X[out_id1, ] %*% Beta # risk score

  # fit our cox model on training data
  fit_cox <- coxph(Surv(time, event) ~ X, data=data_train, x=T)
  
  # in-sample Brier score
  ut_train <- unique(data_train$time[data_train$event==1])
  brier_in <- brier(fit_cox, ut_train, newdata = data_train)
  
  # out of sample Brier score
  ## case 0
  ut_test <- unique(data_test$time[data_test$event==1])
  brier_out0 <- brier(fit_cox, ut_test, newdata = data_test)
 
  ## case 1
  brier_out1 <- brier(fit_cox, ut_test, newdata = data_test1)
 
  ## case 2
  brier_out2 <- brier(fit_cox, ut_test, newdata = data_test2)
  
  # clean results
  tv_bs_df <- bind_rows(
    data.frame(t = ut_train, bs = brier_in$brier, 
               sample = "In-sample", outlier = "None"),
    data.frame(t = ut_test, bs = brier_out0$brier, 
               sample = "Out-of-sample", outlier = "None"),
    data.frame(t = ut_test, bs = brier_out1$brier, 
               sample = "Out-of-sample", outlier = "N(5, 1)"),
    data.frame(t = ut_test, bs = brier_out2$brier, 
               sample = "Out-of-sample", outlier = "N(0, 5)")
  )
  
    # save to final results
    bs_lst[[iter]] <- tv_bs_df
    
    
    setTxtProgressBar(pb, value=iter)
    iter <- iter+1
}

close(pb)
t2 <- Sys.time()

t2-t1

#### 

#### check results ####


bind_rows(bs_lst, .id="iter") %>% 
  # left_join(true_auc_df, by = "time") %>% 
  # pivot_longer(3:5, names_to = "estimator", values_to = "AUC") %>%
  ggplot(aes(x=t, y=bs, col=outlier, linetype = sample))+
  geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
              method = "gam")+
  # geom_line(aes(x = time, y = true_auc), na.rm = T, col = "red", show.legend = F)+
  labs(x="Time", y = "Brier score", col = "Misalignment",
       linetype = "Sample")+
  theme(text = element_text(size = 15), axis.text = element_text(size = 10))+
  # facet_wrap(~estimator)+
  scale_colour_manual(values=cbPalette)
ggsave(filename = here("Code/BrierScore/bs_contam.pdf"),
       width = 6, height = 4, dpi = 300, bg = "white")

#### Results ####
auc_lst2 <- auc_lst
c_lst2 <- c_lst
train_list2 <- train_lst

save(auc_lst2, c_lst2, file = here("Data/SimContamData.RData"))
save(train_list2, file = here("Data/SimContamTrainData.RData"))
