rm(list = ls())
library("survival")
library("ggplot2")
library("cubature")
library("tidyverse")
library(here)
library(risksetROC)
library(mgcv)
library(scam)
# set the seed
set.seed(102131)


# helper functions
source(here("Code/helpers.R"))

#### Simulation set up ####

M <- 1000 # number of simulated data sets
N <- 500 # number of subjects
Beta <- c(1,-1,0.25)

# delta_t <- 0.05

# result container
auc_lst <- list()
c_lst <- list()

##### simulation ####
pb <- txtProgressBar(min=0, max=M,style=3)
for(iter in 1:M){
  # generate data
  X  <- matrix(rnorm(N*3), ncol=3, nrow=N)
  eta <- X %*% Beta
  data <- gen_St(eta=eta,lambda=2, p=2, gen_Ct = function(N) rep(1,N))
  
  # get unique t and eta
  t_vec <- unique(data$time[data$event==1])
  eta_vec <- unique(data$eta[data$event==1])
  nt_pred <- length(t_vec)
  neta_pred <- length(eta_vec)
  
  # results container
  auc_mat <- matrix(NA, nrow = nt_pred, ncol = 4)
  colnames(auc_mat) <- c("time", "HZ", "empirical", "sm_empirical")
  auc_mat[, "time"] <- t_vec    
  
  c_mat <- matrix(NA, nrow = 1, ncol = 7)
  colnames(c_mat) <- c("empirical", "HZ_HZ", "empirical_HZ", "sm_empirical_HZ",
                       "HZ_SmS", "empirical_SmS", "sm_empirical_SmS")
  
  # calculate HZ and empirical time_varying AUC
  for(i in seq_along(t_vec)){
    t_i <- t_vec[i]
    # HZ
    auc_mat[i, "HZ"] <- CoxWeights(marker = data$eta, Stime = data$time, status = data$event,
                                      predict.time = t_i, entry = NULL)$AUC
    # trying to get non-parametric estimator
    auc_mat[i, "empirical"] <- ID_AUC(marker = data$eta, Stime = data$time, status = data$event,
                                  predict.time = t_i,  entry = NULL)
  }
  # smoothed empirical AUC
  sm_fit <- gam(empirical ~ s(time, k = 30), data = as.data.frame(auc_mat))
  auc_mat[, "sm_empirical"] <- predict(sm_fit)
  
  # add to final results
  auc_lst[[iter]] <- auc_mat
  
  # concordance
  ## First, estimate survival probability from Kaplan Meier curve 
  KM_est <- survfit(Surv(time,event)~1, timefix=FALSE,data=data)
  KM_est <- KM_est$surv[KM_est$n.event>0]
  
  auc_sort <-arrange(data.frame(auc_mat), time)
  ## HZ and smS concordance
  c_mat[, "HZ_HZ"] <- intAUC(auc_sort$HZ, auc_sort$time, KM_est, method = "HZ")
  c_mat[, "empirical_HZ"] <- intAUC(auc_sort$empirical, auc_sort$time, KM_est, method = "HZ")
  c_mat[, "sm_empirical_HZ"] <- intAUC(auc_sort$sm_empirical, auc_sort$time, KM_est, method = "HZ")
  c_mat[, "HZ_SmS"] <- intAUC(auc_sort$HZ, auc_sort$time, KM_est, method = "smS")
  c_mat[, "empirical_SmS"] <- intAUC(auc_sort$empirical, auc_sort$time, KM_est, method = "smS")
  c_mat[, "sm_empirical_SmS"] <- intAUC(auc_sort$sm_empirical, auc_sort$time, KM_est, method = "smS")
  c_mat[, "empirical"] <- calc_c(data$eta, data$time, data$event)
  
  c_lst[[iter]] <- c_mat 
  
  setTxtProgressBar(pb, value=iter)
}



# plot 1 iter
auc_lst[[1]] %>% data.frame() %>%
  pivot_longer(cols = 2:4) %>%
  ggplot(aes(y = value, x = as.numeric(time), group = name, col = name))+
  geom_line()

auc_lst <- lapply(auc_lst, as.data.frame)
auc_df <- bind_rows(auc_lst)

## concordance
c_df <- bind_rows(lapply(c_lst, as.data.frame))
c_df <- c_df %>% pivot_longer(1:7, names_to = "estimator", values_to = "concordance")
c_df$estimator <- factor(c_df$estimator, 
                         levels = c("empirical", "HZ_HZ", "HZ_SmS",
                                    "empirical_HZ", "empirical_SmS",
                                    "sm_empirical_HZ", "sm_empirical_SmS"))

save(auc_df, c_df, file = here("outputData/estimated_values.RData"))

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

