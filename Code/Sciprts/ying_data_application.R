## Check for (and install if necessary) packages needed to
##  do data manipulation (tidyverse)
##  do analyses (survey)
##  create tables 1 and 2 (tableone)
pckgs <- c("here","tidyverse")
sapply(pckgs, function(x) if(!require(x,character.only=TRUE,quietly=TRUE)) {
        install.packages(x)
        require(x, character.only=TRUE)
})
rm(list=c("pckgs"))
### Set up data paths
## data path for code, source data creation code
code_path <- here("projects","mortality","code")
## processed_data_path is the main directory where we have our processed accelerometry, covariate, and mortality data
processed_data_path <- here("processed_data")


## read in the mortality data
df_mort <- read_rds(file.path(processed_data_path, "mortality", "Mortality_2015_GH.rds"))
## read in covariate data
df_covar <- read_rds(file.path(processed_data_path, "covariate", "Covariate_GH.rds"))
## merge covariate and mortality data, do some re-coding of factor variables
df_covar <- 
        df_covar %>% 
        left_join(df_mort, by="SEQN") %>% 
        mutate(mortstat = replace(mortstat, ucod_leading %in% "004", 0)) %>% 
        ## replace don't know and refused levels with NA
        ## for all factor variables
        mutate_if(.predicate=is.factor, .funs=function(x) recode_factor(x, "Don't know" = NA_character_, "Refused" = NA_character_)) %>% 
        ## do a bit of extra modifications to alcohol consumption 
        ## and diabetes, create mortality follow-up time in years
        mutate(
                ## create "missing alcohol" factor level 
                alcohol_consumption_fac = factor(alcohol_consumption_fac, levels = c(levels(alcohol_consumption_fac), "Missing alcohol")),
                alcohol_consumption_fac = replace(alcohol_consumption_fac, is.na(alcohol_consumption_fac), "Missing alcohol"),
                ## change borderline diabetics to No diabetes
                diabetes = factor(diabetes, levels=c("No","Yes","Borderline"), labels=c("No","Yes","No")),
                ## create mortality time in years
                event_time_years = permth_exm/12
        ) %>% 
        ## remove individuals who did not participate in the examination portion of the survey
        ## (not eligibile for accelerometry -- should not be included in denominators for calculating
        ##  rates of missing data)
        filter(WTMEC2YR > 0)


## read in subject-level accelerometry features (averaged across days)
df_accel_subj <- read_rds(here(processed_data_path, "accelerometry", "processed_features_subject.rds"))



## merge the subject level data
accel_vars <- c("ASTP","SATP","RA_MIMS","TSleepT",
                "TAT","ABout","SBout","TMVT","TMIMS")
df_analysis_subj <- 
        df_covar %>% 
        full_join(df_accel_subj, by="SEQN") %>% 
        filter(age_years_interview >= 50, age_years_interview < 80) %>% 
        filter(!is.na(BMI), !is.na(PIR)) %>% 
        dplyr::select(event_time_years, mortstat, one_of(paste0(accel_vars,"_mean")),
                      "age_years_interview","BMI","PIR") %>% 
        drop_na()

write_rds(df_analysis_subj,"~/Desktop/Ying_NHANES_application.rds")


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










## function for calculating time-dependent AUC and integrated time-dependent AUC
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
        if("method" %in% names(args)){
                Ct_est <- survfit(Surv(time,1-event)~1, timefix=FALSE,data=data)
                Ct_hat <- unique(Ct_est$surv[Ct_est$n.event > 0])
                
                utimes_c <-unique(Ct_est$time[Ct_est$n.event > 0])
                
                ## get kaplan-meier estimator of the censoring distribution evaluated at each unique survival time
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







library("survival");library("risksetROC");
library(mgcv)

set.seed(18918)
nfolds <- 10
N      <- nrow(df_analysis_subj)
inx_ls <- split(1:N, f=rep(1:nfolds, ceiling(N/nfolds))[1:N])
cv_results_arr <- array(NA, dim=c(3,2,nfolds),
                        dimnames=list("estimator" = c("Heagerty & Zheng",
                                                      "Gonen & Heller",
                                                      "Harrell's C"),
                                      "estimand" = c("In-sample","Out-of-sample"),
                                      "fold"=1:nfolds))
cv_results_arr_lin <- array(NA, dim=c(1,2,nfolds),
                        dimnames=list("estimator" = c("R^2"),
                                      "estimand" = c("In-sample","Out-of-sample"),
                                      "fold"=1:nfolds))
pb <- txtProgressBar(min=0,max=nfolds,style=3)
for(k in 1:nfolds){
        ## separate test and training data
        df_train_k <- df_analysis_subj[-inx_ls[[k]],]
        df_test_k <- df_analysis_subj[inx_ls[[k]],]
        

        
        
        
        fit_lin_k <- try(gam(age_years_interview ~ s(ASTP_mean,RA_MIMS_mean, 
                                                     PIR, BMI, TMIMS_mean, k=200, fx=TRUE), 
                         data=df_train_k))
        fit_k <- try(gam(event_time_years ~ s(ASTP_mean,RA_MIMS_mean, 
                                              age_years_interview, BMI, TMIMS_mean, k=100, fx=TRUE), 
                   family=cox.ph, data=df_train_k, weights=mortstat))
        if(!inherits(fit_lin_k,"try-error")){
                eta_test_lin_k <- predict(fit_lin_k, newdata=df_test_k)
                eta_train_lin_k <- predict(fit_lin_k, newdata=df_train_k)
                
                cv_results_arr_lin[1,1,k] <- 1 - sum((eta_train_lin_k - df_train_k$age_years_interview)^2)/sum(( df_train_k$age_years_interview-mean( df_train_k$age_years_interview))^2)
                cv_results_arr_lin[1,2,k] <- 1 - sum((eta_test_lin_k - df_test_k$age_years_interview)^2)/sum(( df_test_k$age_years_interview-mean( df_test_k$age_years_interview))^2)
        }
        if(!inherits(fit_k,"try-error")){
                eta_train_k <- predict(fit_k, newdata=df_train_k)       
                eta_test_k <- predict(fit_k, newdata=df_test_k)       
                
                ## get KM estimete of survival for the training and test datasets
                KM_fit_train <- survfit(Surv(event_time_years,mortstat) ~ 1, timefix=FALSE, data=df_train_k)
                KM_fit_test  <- survfit(Surv(event_time_years,mortstat) ~ 1, timefix=FALSE, data=df_test_k)
                ## obtain corresponding times, t, and KM estimates of S(t)
                ut_train <- unique(KM_fit_train$time[KM_fit_train$n.event > 0])
                ut_test  <- unique(KM_fit_test$time[KM_fit_test$n.event > 0])
                
                AUC_train <- rep(NA, length(ut_train))
                AUC_test  <- rep(NA, length(ut_test))
                AUC_train_emp <- rep(NA, length(ut_train))
                
                df_emp <- df_train_k
                df_emp$time <- df_train_k$event_time_years
                df_emp$event <- df_train_k$mortstat
                df_emp$eta <- eta_train_k
                df_plt <- get_iAUC(data=df_emp, AUC_fn=ID_AUC, iAUC_fn=intAUC, marker_var=paste0("eta"),returnAUC=TRUE)[[1]]
                for(j in 1:length(ut_train)){
                        AUC_train[j] <- CoxWeights(marker=eta_train_k, Stime=df_train_k$event_time_years, 
                                                   status=df_train_k$mortstat, predict.time=ut_train[j])$AUC
                        
                        # AUC_train_emp[j] <-  get_iAUC(data=df_emp, AUC_fn=intAUC, iAUC_fn=IntegrateAUC, marker_var=paste0("eta"),
                        #                               returnAUC=TRUE)
                        
                        
                }
                for(j in 1:length(ut_test)){
                        AUC_test[j] <- CoxWeights(marker=eta_test_k, Stime=df_test_k$event_time_years, 
                                                  status=df_test_k$mortstat, predict.time=ut_test[j])$AUC
                }
                
                cv_results_arr[1,1,k] <- IntegrateAUC(AUC_train, ut_train, KM_fit_train$surv[KM_fit_train$n.event > 0], tmax=1)
                cv_results_arr[1,2,k] <- IntegrateAUC(AUC_test, ut_test, KM_fit_test$surv[KM_fit_test$n.event > 0], tmax=1)
                
                inx_sub_k <- sample(1:nrow(df_train_k), size=1000, replace=FALSE)
                cv_results_arr[2,1,k] <- calc_c_gh(beta_hat = coef(fit_k), X = predict(fit_k, newdata=df_train_k[inx_sub_k,], type="lpmatrix"))
                cv_results_arr[2,2,k] <- calc_c_gh(beta_hat = coef(fit_k), X = predict(fit_k, newdata=df_test_k, type="lpmatrix"))
                
                cv_results_arr[3,1,k] <- calc_c(marker=eta_train_k, Stime=df_train_k$event_time_years, status=df_train_k$mortstat)
                cv_results_arr[3,2,k] <- calc_c(marker=eta_test_k, Stime=df_test_k$event_time_years, status=df_test_k$mortstat)
                
                
                
        }
        
        setTxtProgressBar(pb, k)
}


df_results_cox <- as.data.frame.table(cv_results_arr)
df_results_lin <- as.data.frame.table(cv_results_arr_lin)

plt_cox <- 
        df_results_cox %>% 
        mutate(NP = ifelse(estimator == "Harrell's C", "Non-parametric","Semi-parametric")) %>% 
        ggplot() + 
        geom_boxplot(aes(x=estimator,y=Freq, fill=NP)) + 
        facet_grid(~estimand) + theme_classic() + xlab("") + 
        ylab("Concordance") + ggtitle("(A) Time-to-Event Outcome: Cox Regression") + 
        theme(legend.position=c(0.9,0.85)) + labs(fill="")

plt_lin <- 
        df_results_lin %>% 
        mutate(estimator="1 - RSS/TSS") %>% 
        ggplot() + 
        geom_boxplot(aes(x=estimator,y=Freq)) + 
        facet_grid(~estimand) + theme_classic() + xlab("") + 
        ylab(expression(R^2)) + ggtitle("(B) Gaussian Outcome: Linear Regression") +
        geom_hline(yintercept=0,col='grey',lty=2)

ggsave(cowplot::plot_grid(plt_cox, plt_lin,ncol=1), file="~/Desktop/ying_presentation_data_application.jpeg",
       height=7.5,width=10)




plt_HZ <- a
df_plt <- 
        df_plt %>% 
        mutate(AUC_HZ = AUC_train) %>% 
        pivot_longer(cols=c("AUC","AUC_HZ")) %>% 
        mutate(estimator = ifelse(name == "AUC", "Non-parametric","Semi-parametric (HZ)")) 
plt_emp <- 
        df_plt %>% 
        ggplot() + 
        geom_point(aes(x=time,y=value)) + 
        facet_grid(~estimator) +
        theme_classic() + ylab(expression({AUC^(I/D)} (t))) + 
        xlab("Time (t)") + 
        geom_smooth(aes(x=time,y=value), method="gam",formula=y~s(x,bs='cr',k=10),se=FALSE)
        


