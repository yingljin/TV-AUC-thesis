rm(list=ls())
library("survival");library("risksetROC"); library("glmnet")
theme_set(theme_minimal())
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set the seed
set.seed(102131)

## Number of subjects used to fit the model
N_obs  <- 300
## Number of subjects used to estimate out-of-sample Concordance
N_test  <- 500
## Total number of subjects to simulate in each iteration (in sample + out of sample)
N_sim_tot  <- N_obs + N_test

## function for simulating survival times with Weibull baseline hazard
# ?independent censor and survival time?
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


## Write a function to calculate Concordance using concordant vs discordant pairs empricical estimator
## Could use survival::concordance, but for some reason I'm getting an error when trying to predict out-of-sample.
calc_c <- function(marker, Stime, status){
    utimes <- sort(unique(Stime[status==1]))
    num <- denom <- 0
    for(ut in seq_along(utimes)){
        ## current time
        ti    <- utimes[ut]
        ## subjects who experienced an event at current time
        inx_i <- which(Stime == ti & status==1)
        ## subjects with observed times beyond event current time
        inx_j <- which(Stime > ti)
        ## number of "cases" and "controls" at time current time
        n_case_t    <- length(inx_i)
        n_control_t <- length(inx_j)
        for(i in seq_along(inx_i)){
            num   <- num + sum( (marker[inx_j] > marker[inx_i[i]] ) ) + 0.5*sum(marker[inx_j] == marker[inx_i[i]])
        }
        denom <- denom + n_case_t*n_control_t
    }
    1-num/denom
}

## calculate concordance under Gonen and heller
calc_c_gh <- function(beta_hat, X){
    N <- nrow(X)
    C <- 0
    for(i in 1:(N-1)){
        for(j in (i+1):N){
            eta_ij <- (X[i,,drop=F] - X[j,,drop=F]) %*% beta_hat
            eta_ji <- (X[j,,drop=F] - X[i,,drop=F]) %*% beta_hat
            if(eta_ij < 0){
                C <- C + 1/(1 + exp(eta_ij))
            } 
            if(eta_ji < 0){
                C <- C + 1/(1+exp(eta_ji))
            }
        }
    }
    2/(N*(N-1))*C
}



## We fit 3 Cox models of the form: log \lambda_ik(t|.) = log \lambda_0(t) + \eta_ik
## for k = 1,2,3 denoting the model
## where the assumed form of \eta_il differs between the 3 models
##
##    Model 1: the true model using 3 predictors: \eta_i1 = X_i1*\beta_1 +  X_i2*\beta_2 +  X_i3*\beta_3
##    Model 2: true model plus 20 predictors which are just noise: \eta_i2 = \eta_i1 + Z_i1 \gamma_1 + ... + Z_{i,20} \gamma_20
##    Model 3: true model plus 100 predictors which are just noise: \eta_i3 = \eta_i1 + Z_i1 \gamma_1 + ... + Z_{i,100} \gamma_100
##
## True beta vector
Beta <- c(1,-1,0.25)



## Simulate the data

## simulate data that are used to obtain true log hazard for each subject
X  <- matrix(rnorm(N_sim_tot*3), ncol=3, nrow=N_sim_tot)
## simulate noise variables for each subject
Z <- matrix(rnorm(N_sim_tot*100), ncol=100, nrow=N_sim_tot)
## obtain true linear predictor (less log hazard) for each subject
eta <- X %*% Beta
## simulate survival times/event indicators for each subject in the test and training data set
data <- gen_St(eta=eta,lambda=2, p=2)
data$X <- I(X)
data$Z <- I(Z)
## separate out the training and test datasets
data_train <- data[1:N_obs,]
data_test  <- data[-c(1:N_obs),]

## fit our 3 models with different number of noise predictors (0, 20, 100)
fit_1 <- coxph(Surv(time, event) ~ X, data=data_train)
fit_2 <- coxph(Surv(time, event) ~ X + Z[,1:20], data=data_train)
fit_3 <- coxph(Surv(time, event) ~ X + Z, data=data_train)
fit_3_L1 <- glmnet(x= model.matrix(~data_train$X + data_train$Z-1), y=Surv(data_train$time, data_train$event), family="cox", alpha=1)
fit_3_L1_cv <-cv.glmnet(x= model.matrix(~data_train$X + data_train$Z-1), y=Surv(data_train$time, data_train$event), family = "cox",alpha=1)
fit_3_L1_coef <- as.vector(coef(fit_3_L1_cv, s="lambda.min"))


## get estimated \hat\eta_ik for k = 1,2,3 models
# for the training data
data_train$eta1 <- fit_1$linear.predictors
data_train$eta2 <- fit_2$linear.predictors
data_train$eta3 <- fit_3$linear.predictors
data_train$eta4 <- model.matrix(~data_train$X + data_train$Z-1) %*% fit_3_L1_coef 
# for the test data
data_test$eta1 <- as.vector(data_test$X %*% coef(fit_1))
data_test$eta2 <- as.vector(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2))
data_test$eta3 <- as.vector(cbind(data_test$X,data_test$Z) %*% coef(fit_3))
data_test$eta4 <- model.matrix(~data_test$X + data_test$Z-1) %*% fit_3_L1_coef


C_ins <- C_oos <- 
    data.frame("Model"= c("True Model","True Model + 20 Noise Vars", "True Model + 100 Noise Vars", "True Model + 100 Noise Vars + L1 Regularization"),
               "C_HZ05" = rep(NA,4),
               "C_Harr" = rep(NA,4),
               "C_GH" = rep(NA, 4))

## get KM estimete of survival for the training and test datasets
KM_fit_train <- survfit(Surv(time,event) ~ 1, timefix=FALSE, data=data_train)
KM_fit_test  <- survfit(Surv(time,event) ~ 1, timefix=FALSE, data=data_test)
## obtain corresponding times, t, and KM estimates of S(t)
ut_train <- unique(KM_fit_train$time[KM_fit_train$n.event > 0])
ut_test  <- unique(KM_fit_test$time[KM_fit_test$n.event > 0])

for(i in 1:4){
    ## create empty vectors for storing AUC(t) 
    AUC_train <- rep(NA, length(ut_train))
    AUC_test  <- rep(NA, length(ut_test))
    
    for(j in 1:length(ut_train)){
        AUC_train[j] <- CoxWeights(marker=data_train[[paste0("eta",i)]], Stime=data_train$time, 
                                   status=data_train$event, predict.time=ut_train[j])$AUC
    }
    for(j in 1:length(ut_test)){
        AUC_test[j] <- CoxWeights(marker=data_test[[paste0("eta",i)]], Stime=data_test$time, 
                                   status=data_test$event, predict.time=ut_test[j])$AUC
    }
    
    
    C_ins$C_HZ05[i] <- IntegrateAUC(AUC_train, ut_train, KM_fit_train$surv[KM_fit_train$n.event > 0], tmax=1)
    C_oos$C_HZ05[i] <- IntegrateAUC(AUC_test, ut_test, KM_fit_test$surv[KM_fit_test$n.event > 0], tmax=1)
    
    C_ins$C_Harr[i] <- calc_c(marker=data_train[[paste0("eta",i)]], Stime=data_train$time, status=data_train$event)
    C_oos$C_Harr[i] <- calc_c(marker=data_test[[paste0("eta",i)]], Stime=data_test$time, status=data_test$event)
    
}

C_ins$C_GH[1] <- calc_c_gh(beta_hat = coef(fit_1), X = unclass(data_test$X))
C_ins$C_GH[2] <- calc_c_gh(beta_hat = coef(fit_2), X = cbind(data_test$X, data_test$Z[,1:20]))
C_ins$C_GH[3] <- calc_c_gh(beta_hat = coef(fit_3), X = cbind(data_test$X, data_test$Z[,1:100]))
C_ins$C_GH[4] <- calc_c_gh(beta_hat = fit_3_L1_coef, X = cbind(data_test$X, data_test$Z[,1:100]))



C_oos$C_GH[1] <- calc_c_gh(beta_hat = coef(fit_1), X = unclass(data_test$X))
C_oos$C_GH[2] <- calc_c_gh(beta_hat = coef(fit_2), X = cbind(data_test$X, data_test$Z[,1:20]))
C_oos$C_GH[3] <- calc_c_gh(beta_hat = coef(fit_3), X = cbind(data_test$X, data_test$Z[,1:100]))
C_oos$C_GH[4] <- calc_c_gh(beta_hat = fit_3_L1_coef, X = cbind(data_test$X, data_test$Z[,1:100]))

## Out-ofpsample C-index using risksetROC increases to 0.90 with 100 noise variables 
C_oos
## In-sample estimates comparable
C_ins













#####################
## Some plots 
####################
library("ggplot2");
library("dplyr")


## To see what's happening, consider fixing t and looking at the ROC curve at a single point 
## when we add a single extreme outlier to the data

# fix time 
ti <- ut_test[50]


# estimate incident/dynamic ROC curve using risksetROC::CoxWeights
# show the effect of outliers on the H&Z estimator by adding a single outlier to the data 
# who lives until toward the end of the study
# remove the first observation
marker_fake <- c(data_train$eta1[-1],  10*sd(data_train$eta1))
time_fake   <- c(data_train$time[-1], 0.85)
event_fake  <- c(data_train$event[-1], 1)
roc_HZ05_eta1 <- CoxWeights(marker=data_train$eta1, Stime=data_train$time, 
                            status=data_train$event, predict.time=ti)
roc_HZ05_eta_fake <- CoxWeights(marker=marker_fake, Stime=time_fake, 
                            status=event_fake, predict.time=ti)

# transform the data to ggplot friendly format
data_plt <- data.frame("Estimator" = paste0("Haegerty and Zheng 2005: t = ", round(ti,2)), 
                          "Model" = c(rep("Original data", length(roc_HZ05_eta1$FP)),
                                      rep("Original data with 1 Outlier", length(roc_HZ05_eta_fake$FP))),
                          "TP" = c(roc_HZ05_eta1$TP,roc_HZ05_eta_fake$TP), 
                          "FP" = c(roc_HZ05_eta1$FP,roc_HZ05_eta_fake$FP))

table(data_plt$Model)

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
ggsave(filename="outlier_exp.png", path="Images/", 
       width=7, height=4, bg="white")

# weight
exp(marker_fake[300])/sum(exp(marker_fake[time_fake>=ti]))


# + 
    ggtitle("Incident/Dynamic ROC Curve at a single time point \n Comparing Estimated Results when Adding in A Single Outlier ")

figure_path <- "~/Desktop"
jpeg(file.path(figure_path, "concordance_1_face_observation.jpeg"),height=400,width=550,quality=100)
    data_plt %>% 
        ggplot() + 
        geom_line(aes(x=FP, y=TP, color=Model)) + 
        facet_grid(. ~ Estimator) +  theme_bw() +
        geom_text(  data=data_text, mapping = aes(x = xind, y =yind, label = label,color=Model)) + 
        xlab("False Positive Rate") + ylab("True Positive Rate") + 
        ggtitle("Incident/Dynamic ROC Curve at a single time point \n Comparing Estimated In-Sample Results when Adding in A Single Fake Outlier ")
dev.off()
