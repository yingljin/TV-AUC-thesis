
##### set-up #####
library("survival");library("risksetROC");
library(mgcv)
library(here)
library(tidyverse)
library(clinfun)
library(scam)

df_analysis_subj <-  read_rds(here("data/Ying_NHANES_application.rds"))


set.seed(18918)
nfolds <- 10

N   <- nrow(df_analysis_subj) # sample size

source(here("Code/Simulation/helpers.R"))
source(here("Code/DataAppl/helpers_appl.R"))

# split fold
inx_ls <- split(1:N, f=rep(1:nfolds, ceiling(N/nfolds))[1:N])

# container for TV-AUC
tvauc_in <- list()
tvauc_out <- list()

# container for concordance
cv_results_arr <- array(NA, dim=c(5,2,nfolds),
                        dimnames=list("estimator" = c("HZ",
                                                      "GH",
                                                      "Harrell",
                                                      "NP",
                                                      "SNP"),
                                      "estimand" = c("In-sample","Out-of-sample"),
                                      "fold"=1:nfolds))

#### simulation ####

pb <- txtProgressBar(min=0,max=nfolds,style=3)
for(k in 1:nfolds){
  
        ## separate test and training data
        df_train_k <- df_analysis_subj[-inx_ls[[k]],]
        df_test_k <- df_analysis_subj[inx_ls[[k]],]
        
        # event time
        t_uni_train <- unique(df_train_k$event_time_years[df_train_k$mortstat==1])
        nt_uni_train <- length(t_uni_train)
        t_uni_test <- unique(df_test_k$event_time_years[df_test_k$mortstat==1])
        nt_uni_test <- length(t_uni_test)
        
        # cox ph survival model
        fit_k <- try(gam(event_time_years ~ s(ASTP_mean,RA_MIMS_mean, 
                                              age_years_interview, BMI, TMIMS_mean, k=100, fx=TRUE), 
                   family=cox.ph, data=df_train_k, weights=mortstat))

  
        if(!inherits(fit_k,"try-error")){
          # TV-AUC
          ## in-sample
          tvauc_in[[k]] <- train_auc(fit_k, 
                                     data = df_train_k %>% 
                                       select(event_time_years, mortstat) %>% 
                                       rename("time" = event_time_years,
                                              "event" = mortstat), 
                                     t = t_uni_train, nt = nt_uni_train)
          ## out-of-sample
          eta_test <- predict(fit_k, newdata = df_test_k, type = "link")
          tvauc_out[[k]] <- test_auc(eta_test, 
                                     data = df_test_k %>% 
                                       select(event_time_years, mortstat) %>% 
                                       rename("time" = event_time_years,
                                              "event" = mortstat),
                                     t = t_uni_test, nt = nt_uni_test)
          
          # concordance
          
          ## get KM estimates of survival for the training and test datasets
          KM_fit_train <- survfit(Surv(event_time_years,mortstat) ~ 1, timefix=FALSE, data=df_train_k)
          KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]
          KM_fit_test  <- survfit(Surv(event_time_years,mortstat) ~ 1, timefix=FALSE, data=df_test_k)
          KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]
          
          ## HZ
          cv_results_arr["HZ","In-sample",k] <- IntegrateAUC(tvauc_in[[k]][, "HZ"], 
                                                             t_uni_train, 
                                                             KM_est_train, tmax=1)
          cv_results_arr["HZ", "Out-of-sample",k] <- IntegrateAUC(tvauc_out[[k]][, "HZ"],
                                                           t_uni_test,
                                                           KM_est_test, tmax=1)
           
         ## Gonan and Heller
         inx_sub_k <- sample(1:nrow(df_train_k), size=1000, replace=FALSE)
         cv_results_arr["GH","In-sample",k] <- calc_c_gh(beta_hat = coef(fit_k),
                                                        X = predict(fit_k, newdata=df_train_k[inx_sub_k,], type="lpmatrix"))
         cv_results_arr["GH","Out-of-sample",k] <- calc_c_gh(beta_hat = coef(fit_k), 
                                                            X = predict(fit_k, newdata=df_test_k, type="lpmatrix"))
        
         ## Harrell's C
         cv_results_arr["Harrell", "In-sample", k] <- calc_c(fit_k$linear.predictors, 
                                                            Stime=df_train_k$event_time_years, 
                                                            status=df_train_k$mortstat)
         cv_results_arr["Harrell", "Out-of-sample", k] <- calc_c(marker=eta_test, 
                                                                Stime=df_test_k$event_time_years,
                                                                status=df_test_k$mortstat)
        ## Non-parametric
        auc_sort_in <-arrange(data.frame(tvauc_in[[k]]), time)
        auc_sort_out <-arrange(data.frame(tvauc_out[[k]]), time)
        cv_results_arr["NP", "In-sample", k] <- intAUC(auc_sort_in$NP, 
                                                       auc_sort_in$time, 
                                                       KM_est_train, 
                                                       method = "smS")
        cv_results_arr["NP", "Out-of-sample", k] <- intAUC_appl(auc_sort_out$NP, 
                                                       auc_sort_out$time, 
                                                       KM_est_test, 
                                                       method = "smS")
        
        ## Smoothed non-parametric
        cv_results_arr["SNP", "In-sample", k] <- intAUC(auc_sort_in$SNP,
                                                        auc_sort_in$time, 
                                                        KM_est_train, 
                                                        method = "smS")
        cv_results_arr["SNP", "Out-of-sample", k] <- intAUC_appl(auc_sort_out$SNP,
                                                        auc_sort_out$time, 
                                                        KM_est_test, 
                                                        method = "smS")
        }
        
        setTxtProgressBar(pb, k)
}

##### results #####

# concordance

df_results_cox <- as.data.frame.table(cv_results_arr)

plt_cox <- df_results_cox %>% 
  mutate(type = ifelse(estimator %in% c("Harrell", "NP","SNP"), 
                       "Non-parametric", "Semi-parametric")) %>% 
        ggplot() + 
        geom_boxplot(aes(x=estimator,y=Freq, fill=type)) + 
        facet_grid(~estimand) + theme_classic() + xlab("") + 
        ylab("Concordance") + ggtitle("(A) Time-to-Event Outcome: Cox Regression") + 
        theme(legend.position=c(0.9,0.85)) + labs(fill="")

# TV-AUC
tvauc_in_df <- lapply(tvauc_in, as.data.frame) %>%
  bind_rows(.id = "iter") %>%
  mutate(sample = "In-sample")

tvauc_out_df <- lapply(tvauc_out, as.data.frame) %>%
  bind_rows(.id = "iter") %>%
  mutate(sample = "Out-of-sample")

bind_rows(tvauc_in_df, tvauc_out_df) %>%
  pivot_longer(3:5, names_to = "Estimator", values_to = "AUC") %>%
  ggplot(aes(x = time, y = AUC, col = sample))+
  geom_smooth(se = F)+
  facet_grid(~Estimator)
  