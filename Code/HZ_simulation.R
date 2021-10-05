##### set up #####
set.seed(925)

library(MASS)
library(tidyverse)
library(survival)
library(mvtnorm)
library(pROC)
library(risksetROC)
library(ggplot2)


##### Generate data #####

# marker values and survival times (logged) are jointly generated from bivariate normal
# rho is correlation
# does that mean covariance is rho^2 ?
gen_df <- function(N = 200, rho = -0.7, mu_binom = c(0, 0), censor_rate = 0.4){
  # biomarker and event time
  sigma_mat = matrix(c(1, -rho, -rho, 1), nrow = 2, ncol = 2)
  df <- mvrnorm(N, mu=mu_binom, Sigma = sigma_mat) %>%
    data.frame() %>%
    rename("M" = X1, "T_event" = X2) %>%
    mutate(T_event = exp(T_event))
  # generate censor time from lognormal distribution
  # first get mean censor time from censor rate
  # log(censor time) ~ normal, log(event time) ~ normal, also independent
  # I've derived this formula myself and check by the censor rate of generated samples
  mu_censor <- mu_binom[2]-qnorm(censor_rate)*sqrt(1+1)
  df$T_censor <- rnorm(N, mu_censor, 1) # censor time generated independently from log-normal distribution
  df$T_censor <- exp(df$T_censor)
  df <- df %>%
    mutate(time = pmin(T_event, T_censor),
           status = as.numeric(T_event < T_censor))
  
  
  # censor rate
  cr <- sum(df$T_event>=df$T_censor)/N
  
  return(list(df = df, cr = cr))
}

# test the function
# df <- gen_df()
# df <- df$df
# df %>% select(M, T_event, T_censor) %>% 
#   mutate_at(vars(T_event, T_censor), log) %>%
#  lapply(summary)
# 
# cov(df %>% select(M, T_event) %>% mutate(T_event = log(T_event)))
# 1-mean(df$status)

#####  calculate true AUC #####

true_auc <- function(df, rho = -0.7, mu_binom = c(0, 0)){
  ut <- unique(df$time[df$status==1]) # unique event time
  c_vec <-  unique(df$M) # unique biomarker values as cutoff sequence
  tauc <- rep(NA, length(ut))
  for(i in seq_along(ut)){
    logt <- log(ut[i])
    # sensitivity: formula from HZ paper in Section 2.5
    true_tp <- pnorm((rho*logt-c_vec)/sqrt(1-rho^2))
    # 1-specificity: formula from HZ paper in Section 2.5
    mat <- matrix(c(c_vec, rep(logt,length(c_vec))), ncol = 2)
    true_fp <- apply(mat, 1, pmvnorm, upper = Inf, mean = mu_binom,
                     corr = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2))/pnorm(-logt)
    # AUC
    id <- order(c_vec)
    true_tp <- true_tp[id]
    true_fp <- true_fp[id]
    dFP <- abs(true_fp[-1] - true_fp[-length(true_fp)])
    aTP <- 0.5 * (true_tp[-1] + true_tp[-length(true_tp)])
    tauc[i] <- sum(dFP * aTP)
  }
  
  return(data.frame(ut = ut, tauc = tauc))
}

# test function
# df_tauc <- true_auc(df)
# 
# plot(log(df_tauc$ut), df_tauc$tauc)

##### get MLE of bivariate nomal distribution of (logT, M) #####
# still working on it

#### calculate AUC based on time varying beta ####
# simply modify risksetROC::CoxWeights
# change eta to beta(t)*marker value
auc_TV_beta <- function(marker, Stime, status, predict.time, beta, entry = NULL){
  
  # estimated marker values using time-varying beta
  Target <- predict.time
  eta <- beta*df$M # calculate estimate eta
  
  if (length(entry) == 0) {
    entry = rep(0, NROW(Stime))
  }
  
  # risk set
  at.risk <- ((Stime >= Target) & (entry <= Target))
  the.eta <- eta[at.risk]
  n <- length(the.eta)
  
  # event set
  the.dead <- (Stime == Target) & (status == 1)
  the.dead <- the.dead[at.risk]
  n.dead <- sum(the.dead)
  
  # weight 
  p0 <- rep(1/(n - n.dead), n)
  if (n.dead > 0) p0[the.dead] <- 0
  p1 <- exp(the.eta)/sum(exp(the.eta)) # proportional hazard assumption? 
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum(dFP * aTP)
  out <- list(marker = eta.out, TP = TP, FP = FP, AUC = area)
  return(out)
}


#### Estimate AUC ####
M <- 10 # number of iteration
N <- 200 # sample size in each iteration
c_rate <- rep(NA, M) # examine censor rate
results <- matrix(NA, ncol = 3)
df_auc <- list() 
# loop
pb = txtProgressBar(min = 0, max = M, initial = 0) 
for(i in 1:M){
  # generate new data
  lst <- gen_df()
  df <- lst$df
  c_rate[i] <- lst$cr
  # save unique event time and true AUC
  results <- true_auc(df)
  results$iter <- i
  ut = unique(df$time[df$status==1]) 
  
  # get time varying beta from local cox
  beta1 <- llCoxReg(df$time, NULL, df$status, df$M)$beta[, 1]
  # get time varying beta from Schoenfeld residual
  cox_fit <- coxph(Surv(time,status)~M, data = df)
  beta2 <- SchoenSmooth(cox_fit, df$time, df$status)$beta
  
  # iterations
  results$mple_cox_auc <- rep(NA, length(ut))
  results$mple_cox_local_auc <- rep(NA, length(ut))
  results$schone_smooth <- rep(NA, length(ut))
  for(j in 1:length(ut)){
    # method 2: maximum partial likelihood using the Cox model
    results$mple_cox_auc[j] <- CoxWeights(marker = df$M, 
                               Stime = df$time, 
                               status = df$status,
                               predict.time = ut[j],
                               entry = NULL)$AUC
    
    # method 3: locally weighted partial likelihood using local cox model
    results$mple_cox_local_auc[j] <- auc_TV_beta(marker = df$M, 
                                                   Stime = df$time, 
                                                   status = df$status,
                                                   predict.time = ut[j],
                                                   beta = beta1[j],
                                                   entry = NULL)$AUC
    
    # method 4: simple local linear smoothing of the scaled Schoenfeld residuals
    results$schone_smooth[j] <- auc_TV_beta(marker = df$M, 
                                            Stime = df$time, 
                                            status = df$status,
                                            predict.time = ut[j],
                                            beta = beta2[j],
                                            entry = NULL)$AUC
    
  }
  
  df_auc[[i]] <- results
  setTxtProgressBar(pb,i)
}
close(pb)


##### visualization #####
df_auc <- bind_rows(df_auc)
df_auc_long <- df_auc %>%
  relocate(iter, .before = 1) %>%
  pivot_longer(cols = 3:6, names_to = "method", values_to = "auc")

# scattered AUC(t)
ggplot(df_auc_long, aes(x = log(ut), y = auc, group = method, col = method))+
  geom_point()

# smoothed AUC(t)
ggplot(df_auc_long, aes(x = log(ut), y = auc, group = method, col = method))+
  geom_smooth(method = "gam", formula = y~s(x, bs = "cr"))
