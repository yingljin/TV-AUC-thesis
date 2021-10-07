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

M <- 100
N <- 200
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
  ## First, estimate survival probablity from Kaplan Meier curve (or cox model?)
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



##### true AUC #####
auc_lst <- lapply(auc_lst, as.data.frame)
auc_df <- bind_rows(auc_lst)
# to speed up, binned time to fewer points
brk <- seq(0, 1, 0.01)
auc_df$time_bin <- cut(auc_df$time, breaks = brk, include.lowest = T,
                       labels = seq(0, 1, length.out = 100))
auc_df$time_bin <- as.numeric(as.character(auc_df$time_bin))

tind <- unique(auc_df$time_bin)
etaind <- seq(-5, 5, length.out = 100)
true_auc <- rep(NA, length(tind))

pb <- txtProgressBar(min=0, max=length(tind), style=3)
for(i in seq_along(tind)){
  t_i <- tind[i]
  #data_i_sens <- subset(data, (time >= t_i) & (time <= (t_i + delta_t)) & (event == 1))
  #data_i_spec <- subset(data, (time >= t_i))
  ## containers
  sens <- rep(NA, length(etaind))
  spec <- rep(NA, length(etaind))
  for(j in seq_along(etaind)){
    eta_ij <- etaind[j]
    
    sens[j] <- adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(eta_ij), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral/
      adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(-Inf), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral
    
    spec[j] <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(eta_ij, Inf), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral/
      adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(Inf, Inf), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral
  }
  ## integrate to AUC
  true_auc[i] <- trap_integrate_ROC(etaind, sens, spec)
  
  setTxtProgressBar(pb, value=i)
}

# brief look at true auc
true_auc <- data.frame(time_bin = tind, estimator = "true", auc = true_auc)

##### true concordance #####
true_auc_sort <- true_auc %>% 
  mutate(time_bin = as.numeric(time_bin)) %>%
  arrange(time_bin)


## marginal S(t)
tind <- true_auc_sort$time_bin
sig_eta <- sum(Beta^2)
true_marg_st <- rep(NA, length(tind))
for(i in seq_along(true_marg_st)){
  true_marg_st[i] <- adaptIntegrate(true_st, t=tind[i], lambda=2, p=2, sigma_eta=sig_eta,
                                    lowerLimit = -5, upperLimit = 5)$integral
}

plot(tind, true_marg_st)

## marginal f(t)
true_marg_ft <- rep(NA, length(tind))
for(i in seq_along(true_marg_ft)){
  true_marg_ft[i] <-adaptIntegrate(true_ft, t=tind[i], lambda=2, p=2, sigma_eta=sig_eta,
                                    lowerLimit = -5, upperLimit = 5)$integral
}

plot(tind, true_marg_ft)


## use trapezoidal rule to approximate integral
y_vec <- 2 * true_marg_ft * true_marg_st*true_auc_sort$auc

plot(tind, y_vec)
  
nt <- nrow(true_auc_sort)
width <- diff(true_auc_sort$time_bin)
height <- y_vec[1:nt-1]+y_vec[2:nt]
true_c <- sum(width*height/2, na.rm = T)


####  plot all iterations #####
auc_df_mean <- auc_df %>% 
  dplyr::select(time_bin, HZ, empirical, sm_empirical) %>%
  pivot_longer(2:4, names_to = "estimator", values_to = "auc") %>%
  group_by(estimator, time_bin) %>% 
  summarize_at("auc", mean)

auc_df_mean <- rbind(auc_df_mean, true_auc)

auc_df_mean %>%
  ggplot(aes(x = time_bin, y = auc, col = estimator))+
  geom_line()+
  labs(title = "Average time-varying AUC")

auc_df[, c(3, 5)] %>%
  ggplot(aes(x = factor(time_bin), y = empirical))+
  geom_boxplot()

auc_df[, c(4, 5)] %>%
  ggplot(aes(x = factor(time_bin), y = sm_empirical))+
  geom_boxplot()
  
#### smooth instead of bin #####

auc_df %>% 
  left_join(true_auc %>% dplyr::select(time_bin, auc) %>% rename(true = auc)) %>%
  relocate(time_bin, .before = 1) %>%
  pivot_longer(3:6, names_to = "estimator", values_to = "auc") %>%
  ggplot(aes(x = time, y = auc, group = estimator, col = estimator))+
  geom_smooth(se = F)
  

##### concordance #####
c_df <- bind_rows(lapply(c_lst, as.data.frame))
c_df <- c_df %>% pivot_longer(1:7, names_to = "estimator", values_to = "concordance")
c_df$estimator <- factor(c_df$estimator, 
                         levels = c("empirical", "HZ_HZ", "HZ_SmS",
                                    "empirical_HZ", "empirical_SmS",
                                    "sm_empirical_HZ", "sm_empirical_SmS"))


c_df %>%  ggplot(aes(x = estimator, y = concordance))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 60))+
  geom_hline(yintercept = true_c)


