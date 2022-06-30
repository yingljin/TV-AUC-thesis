rm(list=ls())
library("risksetROC"); 
library("mgcv");       # package used for smoothing AUC(t)
library("scam");       # package used for fitting a monotonically increasing penalized spline model to the marginal survival curve

## Number of subjects used to fit the model
N_obs  <- 300
## Number of subjects used to estimate out-of-sample Concordance
N_test     <- c(500, 1500)
## Total number of subjects to simulate in each iteration (in sample + out of sample)
N_sim_tot  <- N_obs + sum(N_test)
## number of datasets to simulate
nsim <- 100

## function for simulating survival times with Weibull baseline hazard
# eta: vector containing the log hazard (linear predictor) for each subject
# lambda: scale parameter for the baseline hazard
# p: shape paramter for the baseline hazard
# gen_Ct: function for generating censoring times, default is exponential with mean 5 truncated at t=1
gen_St <- function(eta, lambda, p, gen_Ct = function(N) pmin(1,rexp(N, 1/5)) ){
        N <- length(eta)
        U <- runif(N, 0, 1)
        St = (-log(U)/(lambda*exp(eta)))^(1/p)
        Ct = gen_Ct(N)
        
        data.frame("event" = as.numeric(St<=Ct), "time" = pmin(St,Ct), "eta" = eta)
}



#### function for calculating time-dependent AUC and integrated time-dependent AUC ####
# data: data frame with all information needed to calculate time dependent AUC
# AUC_fn: function which calculates AUC(t)
# iAUC_fn: function which calculates concordance
# returnAUC: set to TRUE if you want both AUC(t) and Concordance returned
get_iAUC <- function(data, AUC_fn, iAUC_fn, returnAUC=FALSE,
                     marker_var="eta", time_var="time", event_var="event",
                     smoothAUC=FALSE,...){
        args <- c(...)
        Ct = NULL
        n_event = NULL
        utimes <- unique(sort(data$time[data$event==1]))
        KM_est <- survfit(Surv(time,event)~1, timefix=FALSE,data=data)
        KM_est <- unique(KM_est$surv[KM_est$n.event > 0])
        
        # question: what is this chunk actually doing?
        if("method" %in% names(args)){
            Ct_est <- survfit(Surv(time,1-event)~1, timefix=FALSE,data=data)
            Ct_hat <- unique(Ct_est$surv[Ct_est$n.event > 0])
            
            utimes_c <-unique(Ct_est$time[Ct_est$n.event > 0])
            
            ## get kaplan-meier estimator of the censoring distribution evaluated at each unique survival time
            ## utimes: unique event time
            ## utimes_c: unique censor time
            Ct   <- vapply(utimes, function(x){
                    inx_x <- which(utimes_c < x)
                    if(length(inx_x) == 0){
                            return(1)
                    } else {
                            return(Ct_hat[max(inx_x)])
                    }
                }, numeric(1))
            
            
            n_event <- vapply(utimes, function(x) sum(data$time == x * data$event), numeric(1))

        }
        # question end
        n_ut    <- length(utimes)
        AUC_vec <- rep(NA, n_ut)
        for(j in seq_along(1:n_ut)){
                AUC_vec[j] <- AUC_fn(marker = data[[marker_var]], 
                                     Stime=data[[time_var]], 
                                     status=data[[event_var]],
                                     predict.time=utimes[j])
        }
        if(smoothAUC){
            fit_AUC <- gam(AUC_vec ~ s(utimes, k=30, bs="tp"))
            AUC_vec <- fit_AUC$fitted.values
        }
        if(is.null(Ct)){
            iAUC <- iAUC_fn(AUC=AUC_vec, utimes=utimes, St=KM_est, tmax=Inf,...)
        } else {
            iAUC <- iAUC_fn(AUC=AUC_vec, utimes=utimes, St=KM_est, tmax=Inf, Ct=Ct, n_event=n_event,...)    
        }
    
        
        if(returnAUC){
                out <- list("AUC" = data.frame("time" = utimes,"AUC" = AUC_vec),
                            "iAUC" = iAUC)
        } else {
                out <- list("iAUC" = iAUC)
        }
        
        out
}


## a function which modifies risksetROC::CoxWeights 
## to estimate time-dependent AUC using an empirical estimator 
ID_AUC <- function(marker, Stime, status, predict.time, entry = NULL, ...){
        if (length(entry) == 0) {
                entry = rep(0, NROW(Stime))
        }
        at.risk <- ((Stime >= predict.time) & (entry <= predict.time))
        eta     <- marker[at.risk]
        status2 <- status
        status2[Stime > predict.time] <- 0
        status2 <- status2[at.risk]
        
        C_num <- 0
        n_case_t    <- sum(status2)
        n_control_t <- sum(1-status2)
        inx_ti <- which(status2 == 1)
        inx_tj <- which(status2 == 0)
        for(id in inx_ti){
                C_num <- C_num + sum(eta[id] > eta[inx_tj]) + 0.5*sum(eta[id] == eta[inx_tj])
        }
        out <- C_num/(n_case_t*n_control_t)
        
        out     
}


## integrate I/D AUC to get concordance measure (C-index)
## IntegrateAUC() in risksetROC assumes an equally-spaced survival time, which may not hold
## So here we create a new function, intAUC(), which uses the trapezoidal rule to peform numeric integration
## Arguments:
##     AUC - vector of I/D auc at each time index in the vector `ut`. 
##     ut - vector of time indices where I/D `AUC` is measured
##     time - vector of event times for all participants
##     event - vector of event indicators for all participants
##     St - vector of marginal survival probabilities at each time index in the vecto `ut` (i.e. not conditioned on any covariates)
##     method - One of the following: HZ, smCH, smS, MC. Note that HZ, smCH, and smS all approximate f(t) under certain assumptions and then numerically 
##              integrate using the trapezoidal rule. The MC option estimates the integral using a monte carlo approximation. 
##              We can make smCH and smS more rigorous by actually calculating the derivative of the estimated curve AND
##              fitting a binomial model (differing number of participants in the risk set at each time)
##        * HZ: this is the method of Haegerty and Zhang. Uses a linear approximation for f(t) 
##        * smS: fits a monotonically increasing B-spline basis to logit of pr(death). We then back out the marginal survival function.
##     smoothAUC - whether to smooth AUC(t) using penalized splines before integrating
intAUC <- function(AUC, utimes, St, method="HZ", smoothAUC=FALSE, n_event=NULL, Ct=NULL, k=30,...){
    ut <- utimes
    if(method %in% c("HZ","smS")){
        if(method == "HZ"){
            ft <- rep(NA, length(St))
            ft[1] <- (1 - St[1])/(ut[1])
            for(j in 2:length(St)){
                ft[j] <- (St[j - 1] - St[j])/(ut[j]-ut[j-1])
            }
        }
        if(method == "smS"){
            fit_S <- scam(St ~ s(ut, bs="mpd",k=k),
                                data=data.frame(ut=ut, St=St))
            df_pred <- data.frame(ut = rep(ut, each=2) + rep(c(-0.00001,0.00001), length(ut)))
            St_pred <- predict(fit_S, newdata=df_pred, type='response')
            ft     <- -diff(St_pred)[seq(1,2*length(ut),by=2)]                    
        }
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
    }
    
    if(method=="MC"){
        if(any(is.null(n_event), is.null(Ct))) stop("If method == \"MC\" must supply n_event and Ct")
        iAUC <- sum(n_event*AUC*St/Ct)/sum(n_event*St/Ct)   
    }
    
    
    return(iAUC)
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





## We fit 3 Cox models of the form: log \lambda_ik(t|.) = log \lambda_0(t) + \eta_ik
## for k = 1,2,3 denoing the model
## where the assumed form of \eta_il differs between the 3 models
##
##    Model 1: the true model using 3 predictors: \eta_i1 = X_i1*\beta_1 +  X_i2*\beta_2 +  X_i3*\beta_3
##    Model 2: true model plus 20 predictors which are just noise: \eta_i2 = \eta_i1 + Z_i1 \gamma_1 + ... + Z_{i,20} \gamma_20
##    Model 3: true model plus 100 predictors which are just noise: \eta_i3 = \eta_i1 + Z_i1 \gamma_1 + ... + Z_{i,100} \gamma_100
##
## True beta vector
Beta <- c(1,-1,0.25)


## We then obtain in-sample and out-of-sample C-index estimates for each model fit using different number of 
## noise predictors, and store them in an matrices. Syntax:
##   - In-sample vs out of sample:
##         - C_*_ins = in-sample 
##         - C_*_oos = out-of-sample)
##   - Different estimators
##         - C_HZ_* = concordance estimated via integrated AUC(t) estimated by risksetROC package 
##         - C_HZ2_* = concordance estimated via integrated AUC(t) estimated by a modified version of the risksetROC::IntegrateAUC function 
##         - C_AUC_* = concordance estimated via integrated AUC(t) estimated using empirical estimator 
##         - C_AUC_* = concordance estimated via integrated AUC(t) estimated using smoothed empirical estimator 
##         - C_MC_* = concordance estimated via integrated AUC(t) estimated Monte-Carlo approach 
##         - C_pairs_* = concordance estimated using the concordant vs discordant pairs estimator
## The rows of the matrix corresponds to simulation iteration
## The columns of the matrix corresponds to different estimators 
##         - [,1] = "True model" (no noise predictors)
##         - [,2] = True model + 20 noise predictors
##         - [,3] = True model + 100 noise predictors
C_HZ_ins <- C_HZ2_ins <- C_AUC_ins <- C_AUC_sm_ins <- C_MC_ins  <- C_pairs_ins <- matrix(NA, nrow=nsim, ncol=3)
C_HZ_oos <- C_HZ2_oos <- C_AUC_oos <- C_AUC_sm_oos <- C_MC_oos <- C_pairs_oos <- array(NA, c(nsim, 3, 2))

## Do the simulation
pb <- txtProgressBar(min=0,max=nsim, style=3)
set.seed(102131)
for(i in 1:nsim){
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
        
        ## get estimated \hat\eta_ik for k = 1,2,3 models
        # for the training data
        data_train$eta1 <- fit_1$linear.predictors
        data_train$eta2 <- fit_2$linear.predictors
        data_train$eta3 <- fit_3$linear.predictors
        # for the test data
        data_test$eta1 <- as.vector(data_test$X %*% coef(fit_1))
        data_test$eta2 <- as.vector(cbind(data_test$X,data_test$Z[,1:20]) %*% coef(fit_2))
        data_test$eta3 <- as.vector(cbind(data_test$X,data_test$Z) %*% coef(fit_3))

        
        for(j in 1:3){
            C_HZ_ins[i,j]    <- get_iAUC(data=data_train, AUC_fn=function(...) CoxWeights(...)$AUC, iAUC_fn=IntegrateAUC, marker_var=paste0("eta",j))$iAUC    
            C_HZ2_ins[i,j]   <- get_iAUC(data=data_train, AUC_fn=function(...) CoxWeights(...)$AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j),method="HZ")$iAUC    
            C_AUC_ins[i,j]   <- get_iAUC(data=data_train,  AUC_fn=ID_AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j), method="smS")$iAUC
            C_AUC_sm_ins[i,j]<- get_iAUC(data=data_train,  AUC_fn=ID_AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j),method="smS",smoothAUC=TRUE)$iAUC
            C_MC_ins[i,j]    <- get_iAUC(data=data_train,  AUC_fn=ID_AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j), method="MC")$iAUC
            C_pairs_ins[i,j] <- calc_c(marker=data_train[[paste0("eta",j)]], Stime=data_train$time, status=data_train$event)
            for(k in 1:2){
                if(k == 1){
                    inx_k <- 1:N_test[1]
                } else {
                    inx_k <- setdiff(1:sum(N_test), 1:N_test[k-1])
                }
                C_HZ_oos[i,j,k]    <- get_iAUC(data=data_test[inx_k,], AUC_fn=function(...) CoxWeights(...)$AUC, iAUC_fn=IntegrateAUC, marker_var=paste0("eta",j))$iAUC
                C_HZ2_oos[i,j,k]   <- get_iAUC(data=data_test[inx_k,], AUC_fn=function(...) CoxWeights(...)$AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j), method="HZ")$iAUC
                C_AUC_oos[i,j,k]   <- get_iAUC(data=data_test[inx_k,],  AUC_fn=ID_AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j),method="smS" )$iAUC
                C_AUC_sm_oos[i,j,k]<- get_iAUC(data=data_test[inx_k,],  AUC_fn=ID_AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j),method="smS",smoothAUC=TRUE)$iAUC
                C_MC_oos[i,j,k]    <- get_iAUC(data=data_test[inx_k,],  AUC_fn=ID_AUC, iAUC_fn=intAUC, marker_var=paste0("eta",j),method="MC" )$iAUC
                C_pairs_oos[i,j,k] <- calc_c(marker=data_test[[paste0("eta",j)]][inx_k], Stime=data_test$time[inx_k], status=data_test$event[inx_k])
            }
        }
        
        setTxtProgressBar(pb, i)
}


library("dplyr");library("ggplot2");library("tidyr"); library("stringr")

results <- data.frame(iteration = rep(rep(1:nsim,3), 2),
                      Model = rep(rep(c(1:3), each=nsim), 2), 
                      N_oos = rep(N_test, each=nsim*3), 
                      C_HZ_ins = rep(as.vector(C_HZ_ins),2),
                      C_HZ2_ins = rep(as.vector(C_HZ2_ins),2),
                      C_AUC_ins = rep(as.vector(C_AUC_ins),2),
                      C_AUC_sm_ins = rep(as.vector(C_AUC_sm_ins),2),
                      C_MC_ins = rep(as.vector(C_MC_ins),2),
                      C_pairs_ins = rep(as.vector(C_pairs_ins),2),
                      C_HZ_oos = as.vector(C_HZ_oos),
                      C_HZ2_oos = as.vector(C_HZ2_oos),
                      C_AUC_oos = as.vector(C_AUC_oos),
                      C_AUC_sm_oos = as.vector(C_AUC_sm_oos),
                      C_MC_oos = as.vector(C_MC_oos),
                      C_pairs_oos =as.vector(C_pairs_oos))

results %>% 
    group_by(Model, N_oos) %>% 
    summarise_all(.funs=mean, na.rm=TRUE)

results <- 
    results %>% 
    gather(Estimator, Concordance, -iteration, -N_oos, -Model) %>% 
    mutate(InSample = as.numeric(grepl("oos", Estimator)),
           Estimator = gsub("\\_oos","",Estimator),
           Estimator = gsub("\\_ins","",Estimator),
           Estimator = gsub("C\\_","",Estimator),
           Estimator = factor(Estimator,
                              levels=c("HZ","HZ2","MC","AUC","AUsm","pairs"),
                              labels=c("Haegerty + Zhang (2005)",
                                       "Haegerty + Zhang -- Corrected(2005)",
                                       "Monte Carlo",
                                       "Integrate Empirical AUC(t)",
                                       "Integrate Smoothed Empirical AUC(t)",
                                       "Empirical Concordance")),
                              # levels=c("HZ","AUC","pairs"),
                              # labels=c("Haegerty + Zhang (2005)",
                              #          "Integrate Empirical AUC(t)",
                              #          "Empirical Concordance")),
           InSample = factor(InSample, levels=c(1,0), labels=c("Out-of-Sample","In-Sample")),
           Model = factor(Model, levels=1:3, labels=c("True Model","20 Noise Variables","100 Noise Variables"))
           )





figure_path <- "~/Desktop"


jpeg(file.path(figure_path, "in_sample_concordance.jpeg",height=1000,width=1100,quality=100))
results %>% 
    ggplot() + 
    geom_boxplot(aes(x=Model, y=Concordance, fill=Estimator)) + 
    facet_grid(InSample~N_oos)
dev.off()

