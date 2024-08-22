
# This script implements the simulation study 
# investigating the effect of overfit
# corresponding to Section 4.1 of the manuscript

##### set up #####
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
set.seed(510)

# helper functions
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

##### simulation ####

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
  
  fit_3 <- tryCatch(expr = coxph(Surv(time, event) ~ X + Z, data=data_train),
                    warning = function(x){NULL})
  
  # if all of the models covereged
  if(!is_null(fit_1) & !is.null(fit_2) & !is.null(fit_3)){
  
    # in-sample BS
    ut_train <- unique(data_train$time[data_train$event==1])
    brier_in1 <- brier(fit_1, ut_train, newdata = data_train)
    brier_in2 <- brier(fit_2, ut_train, newdata = data_train)
    brier_in3 <- brier(fit_3, ut_train, newdata = data_train)
  
    # out-of-sample BS
    ## estimated risk scores
    ut_test <- unique(data_test$time[data_test$event==1])
    brier_out1 <- brier(fit_1, ut_test, newdata = data_test)
    brier_out2 <- brier(fit_2, ut_test, newdata = data_test)
    brier_out3 <- brier(fit_3, ut_test, newdata = data_test)

    
    # put in and out-of-sample together
    tv_bs_df <- bind_rows(
      data.frame(t = ut_train, 
                 bs1 = brier_in1$brier, 
                 bs2 = brier_in2$brier, 
                 bs3 = brier_in3$brier, 
                 sample = "In-sample"),
      
      data.frame(t = ut_test, 
                 bs1 = brier_out1$brier, 
                 bs2 = brier_out2$brier, 
                 bs3 = brier_out3$brier, 
                 sample = "Out-of-sample")
    )
    
    
    # save to final results
    bs_lst[[iter]] <- tv_bs_df
    
    
    # move to next iter
    iter <- iter + 1
  
  }
  # skip non-converged data
  else{skip <- skip+1}

  setTxtProgressBar(pb, value=iter)
}


close(pb)
t2 <- Sys.time()

##### Check results #####



# trend of TV-AUC estimates (smoothed)
bind_rows(bs_lst, .id="iter") %>%
# tv_bs_df %>% 
  # left_join(true_auc_df, by = "time") %>% 
  pivot_longer(3:5, names_to = "model", values_to = "bs") %>%
  # mutate(model = as.factor(model)) %>%
  mutate(model = factor(model, levels = c("bs1", "bs2", "bs3"),
                        labels = c("No noise", "20 noise", "100 noise"))) %>%
  ggplot(aes(x=t, y=bs, col=model, linetype = sample))+
  geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
              method = "gam")+
  # geom_line(aes(x = time, y = true_auc), na.rm = T, col = "red", show.legend = F)+
  labs(x="Time", y = "Brier score", col = "Model", linetype = "Sample")+
  theme(text = element_text(size = 15), axis.text = element_text(size = 10))+
  scale_colour_manual(values=cbPalette)
ggsave(filename = here("Code/BrierScore/bs_overfit.pdf"),
       width = 6, height = 4, dpi = 300, bg = "white")

# concordance
c_df %>% pivot_longer(2:6, names_to = "estimator", values_to = "c") %>%
  ggplot()+
  geom_boxplot(aes(x=sample, y=c))+
  facet_grid(rows=vars(estimator), cols=vars(model))

#### check the NA AUCs ####



for(i in seq_along(auc_lst)){
  
  if(sum(is.na(auc_lst[[i]]$NP)) > 0){
    print(i)
  }
  
}

auc_lst[[4]] %>% filter(is.na(NP))

train_list[[4]] %>% filter(time==0.9708767)
View(train_list[[4]] %>% filter(event == 1) %>% arrange(desc(time)))

#### Save data ####

# save(auc_df_train, auc_df_test, c_df_train, c_df_test, file = here("outputData/estimated_values.RData"))
save(auc_lst, c_lst, file = here("Data/SimNoiseModel.RData"))
save(train_list, file = here("Data/SimNoiseTrainData.RData"))
# load(here("Data/results_by_iter.RData"))

gamma(1+0.5)/sqrt(2)
