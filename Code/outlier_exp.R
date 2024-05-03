
# This script produces the example of one outlier
# corresponding to Figure 1 in the manuscript

library("survival")
library("risksetROC")
library("glmnet")
library(ggplot2)
library("dplyr")
theme_set(theme_minimal())
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set the seed
set.seed(102131)



#### Function for simulating survival times ####
# Weibull baseline hazard
# eta: vector containing the log hazard (linear predictor) for each subject
# lambda: scale parameter for the baseline hazard
# p: shape parameter for the baseline hazard
# gen_Ct: function for generating censoring times, default is exponential with mean 5 truncated at t=1

gen_St <- function(eta, lambda, p, gen_Ct = function(N) pmin(1,rexp(N, 1/5)) ){
    N <- length(eta)
    U <- runif(N, 0, 1)
    St = (-log(U)/(lambda*exp(eta)))^(1/p)
    Ct = gen_Ct(N)
    
    data.frame("event" = as.numeric(St<=Ct), "time" = pmin(St,Ct), "eta" = eta)
}


#### Simulate the data ###

## Number of subjects used to fit the model
N_obs  <- 300
## Number of subjects used for out-of-sample evaluation
N_test  <- 500
## Total number of subjects to simulate 
N_sim_tot  <- N_obs + N_test

## True beta vector
Beta <- c(1,-1,0.25)

## simulate data that are used to obtain true log hazard for each subject
X  <- matrix(rnorm(N_sim_tot*3), ncol=3, nrow=N_sim_tot)

## obtain true linear predictor (less log hazard) for each subject
eta <- X %*% Beta

## simulate survival times/event indicators for each subject in the test and training data set
data <- gen_St(eta=eta,lambda=2, p=2)
data$X <- I(X)

## separate out the training and test datasets
data_train <- data[1:N_obs,]
data_test  <- data[-c(1:N_obs),]

#### Model fit #####

## fit model on training data
fit_1 <- coxph(Surv(time, event) ~ X, data=data_train)

## estimate risk score/linear predictor for the test data
data_test$eta1 <- as.vector(data_test$X %*% coef(fit_1))


#### Out-of-sample AUC ####

# fix time 
ut_test <- unique(data_test$time[data_test$event==1])
ti <- ut_test[50]

# duplicate the data set and introduce one outlier
data_test2 <- data_test
data_test2$X[500, ] <- 40*apply(data_test$X, 2, sd)
data_test2$eta1 <- as.vector(data_test2$X %*% coef(fit_1))
data_test2$time[500] <- 0.85
data_test2$event[500] <- 0

# subjects at risk at t=0.27
summary(survfit(Surv(time, event) ~ 1, data=data_test2))
survfit(fit_1)

data_test[500,]

# survival probablity
est_log_hz <- predict(fit_1, newdata = data_test2)
est_log_hz <- sort(est_log_hz)
est_log_hz[500]/est_log_hz[499]

# calculate AUC
roc_HZ05_eta1 <- CoxWeights(marker=data_test$eta1, Stime=data_test$time, 
                            status=data_test$event, predict.time=ti)
roc_HZ05_eta_fake <- CoxWeights(marker=data_test2$eta1, Stime=data_test2$time, 
                            status=data_test2$event, predict.time=ti)

#### Plot #####
# transform the data to ggplot friendly format
data_plt <- data.frame("Estimator" = paste0("Haegerty and Zheng 2005: t = ", round(ti,2)), 
                          "Model" = c(rep("Original data", length(roc_HZ05_eta1$FP)),
                                      rep("Original data with 1 Outlier", length(roc_HZ05_eta_fake$FP))),
                          "TP" = c(roc_HZ05_eta1$TP,roc_HZ05_eta_fake$TP), 
                          "FP" = c(roc_HZ05_eta1$FP,roc_HZ05_eta_fake$FP))


# separate data frame for adding text to the plot (estimated AUC)
data_text <- data.frame("Estimator" = rep(paste0("Haegerty and Zheng 2005: t = ", round(ti,2)),2),
                        "Model" = c("Original data", "Original data with 1 Outlier"), 
                        "label" = paste0("AUC = ",round(c(roc_HZ05_eta1$AUC, roc_HZ05_eta_fake$AUC),3)),
                        "xind" = rep(0.75, 2),
                        "yind" = c(0.15,0.25))


data_plt %>% 
    ggplot() + 
    geom_line(aes(x=FP, y=TP, color=Model)) + 
    theme_bw() +
    geom_text(data=data_text, mapping = aes(x = xind, y =yind, label = label,color=Model)) + 
    xlab("False Positive Rate") + ylab("True Positive Rate")+
    labs(col="Data")+
  scale_color_manual(values = cbPalette)+
  theme(text=element_text(size=15),
        axis.text = element_text(size=10))
ggsave(filename="Images/outlier_exp.pdf",
       width=7, height=4, bg="white", dpi = 300)
