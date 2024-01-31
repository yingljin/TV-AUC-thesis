
# This code reproduces the data application section 
# in the manuscript

##### set-up #####
rm(list = ls())
library("survival");library("risksetROC");
library(mgcv)
library(here)
library(tidyverse)
library(clinfun)
library(scam)

theme_set(theme_minimal())
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


set.seed(825)


##### data ######

# new updated outcome
df_analysis_subj <-  read_rds(here("Data/NHANES_application_0824.rds"))

# select variables of interest and center & scale covariates
df_analysis_subj <- df_analysis_subj %>% 
  select(event_time_years, mortstat, ASTP_mean, RA_MIMS_mean, age_years_interview,
         BMI, TMIMS_mean) %>%
  mutate_at(vars(ASTP_mean, RA_MIMS_mean, age_years_interview,
                 BMI, TMIMS_mean), scale)

##### source code #####

# same estimators from the simulation section
source(here("Code/Simulation/helpers.R")) 
# functions for saving/formatting results, and calculation of GH 
source(here("Code/DataAppl/helpers_appl.R"))


#####  split fold ####
nfolds <- 10
N   <- nrow(df_analysis_subj) # sample size
inx_ls <- split(1:N, f=rep(1:nfolds, ceiling(N/nfolds))[1:N])


#####  container for TV-AUC #####
# gam-cox
tvauc_in_gam <- list()
tvauc_out_gam <- list()
# linear-cox
tvauc_in_lin <- list()
tvauc_out_lin <- list()

##### container for concordance #####
# gam-cox
c_gam <- array(NA, dim=c(5,2,nfolds),
                        dimnames=list("estimator" = c("HZ",
                                                      "GH",
                                                      "Harrell",
                                                      "NP",
                                                      "SNP"),
                                      "estimand" = c("In-sample","Out-of-sample"),
                                      "fold"=1:nfolds))
# linear-cox
c_lin <- array(NA, dim=c(5,2,nfolds),
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
        
        ## event time
        t_uni_train <- unique(df_train_k$event_time_years[df_train_k$mortstat==1])
        nt_uni_train <- length(t_uni_train)
        t_uni_test <- unique(df_test_k$event_time_years[df_test_k$mortstat==1])
        nt_uni_test <- length(t_uni_test)
        
        # training models
        ## cox-gam
        fit_k_gam <- try(gam(event_time_years ~ s(ASTP_mean,RA_MIMS_mean, 
                                              age_years_interview, BMI, TMIMS_mean, k=200, fx=TRUE), 
                   family=cox.ph, data=df_train_k, weights=mortstat))
        ## cox_linea
        fit_k_lin <- coxph(Surv(event_time_years, mortstat)~ 
                             ASTP_mean + RA_MIMS_mean + age_years_interview + BMI + TMIMS_mean,
                           data = df_train_k)

        # check if there is fitting issues with GAM model
        if(!inherits(fit_k_gam,"try-error")){
          # TV-AUC
          ## in-sample
          tvauc_in_gam[[k]] <- train_auc(fit_k_gam, 
                                     data = df_train_k %>% 
                                       select(event_time_years, mortstat) %>% 
                                       rename("time" = event_time_years,
                                              "event" = mortstat), 
                                     t = t_uni_train, nt = nt_uni_train)
          tvauc_in_lin[[k]] <- train_auc(fit_k_lin, 
                                         data = df_train_k %>% 
                                           select(event_time_years, mortstat) %>% 
                                           rename("time" = event_time_years,
                                                  "event" = mortstat), 
                                         t = t_uni_train, nt = nt_uni_train)
          
          
          ## out-of-sample
          eta_test_gam <- predict(fit_k_gam, newdata = df_test_k, type = "link")
          tvauc_out_gam[[k]] <- test_auc(eta_test_gam, 
                                     data = df_test_k %>% 
                                       dplyr::select(event_time_years, mortstat) %>% 
                                       rename("time" = event_time_years,
                                              "event" = mortstat),
                                     t = t_uni_test, nt = nt_uni_test)
          eta_test_lin <- predict(fit_k_lin, newdata = df_test_k, type = "lp")
          tvauc_out_lin[[k]] <- test_auc(eta_test_lin, 
                                         data = df_test_k %>% 
                                           dplyr::select(event_time_years, mortstat) %>% 
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
          c_gam["HZ","In-sample",k] <- IntegrateAUC(tvauc_in_gam[[k]][, "HZ"], 
                                                    t_uni_train, KM_est_train, tmax=1)
          c_gam["HZ", "Out-of-sample",k] <- IntegrateAUC(tvauc_out_gam[[k]][, "HZ"],
                                                           t_uni_test, KM_est_test, tmax=1)
          c_lin["HZ","In-sample",k] <- IntegrateAUC(tvauc_in_lin[[k]][, "HZ"], 
                                                    t_uni_train, KM_est_train, tmax=1)
          c_lin["HZ", "Out-of-sample",k] <- IntegrateAUC(tvauc_out_lin[[k]][, "HZ"],
                                                         t_uni_test, KM_est_test, tmax=1)
           
         ## Gonan and Heller
         c_gam["GH","In-sample",k] <- calc_c_gh(beta_hat = coef(fit_k_gam),
                                                X = predict(fit_k_gam, newdata=df_train_k, type="lpmatrix"))
         c_gam["GH","Out-of-sample",k] <- calc_c_gh(beta_hat = coef(fit_k_gam), 
                                                    X = predict(fit_k_gam, newdata=df_test_k, type="lpmatrix"))
         c_lin["GH","In-sample",k] <- calc_c_gh(beta_hat = coef(fit_k_lin),
                                                X = df_train_k %>% 
                                                  select(-event_time_years, -mortstat) %>% as.matrix())
         c_lin["GH","Out-of-sample",k] <- calc_c_gh(beta_hat = coef(fit_k_lin), 
                                                    X = df_test_k %>% 
                                                      select(-event_time_years, -mortstat) %>% as.matrix())
         
         ## Harrell's C
         c_gam["Harrell", "In-sample", k] <- calc_c(fit_k_gam$linear.predictors, 
                                                    Stime=df_train_k$event_time_years, 
                                                    status=df_train_k$mortstat)
         c_gam["Harrell", "Out-of-sample", k] <- calc_c(marker=eta_test_gam, 
                                                                Stime=df_test_k$event_time_years,
                                                                status=df_test_k$mortstat)
         c_lin["Harrell", "In-sample", k] <- calc_c(fit_k_lin$linear.predictors,
                                                    Stime=df_train_k$event_time_years, 
                                                    status=df_train_k$mortstat)
         c_lin["Harrell", "Out-of-sample", k] <- calc_c(marker=eta_test_lin, 
                                                        Stime=df_test_k$event_time_years,
                                                        status=df_test_k$mortstat)
        ## Non-parametric
        auc_sort_in_gam <-arrange(data.frame(tvauc_in_gam[[k]]), time)
        auc_sort_out_gam <-arrange(data.frame(tvauc_out_gam[[k]]), time)
        auc_sort_in_lin <-arrange(data.frame(tvauc_in_lin[[k]]), time)
        auc_sort_out_lin <-arrange(data.frame(tvauc_out_lin[[k]]), time)
        c_gam["NP", "In-sample", k] <- intAUC(auc_sort_in_gam$NP, 
                                              auc_sort_in_gam$time, 
                                              KM_est_train, method = "smS")
        c_gam["NP", "Out-of-sample", k] <- intAUC_appl(auc_sort_out_gam$NP, 
                                                       auc_sort_out_gam$time, 
                                                       KM_est_test, 
                                                       method = "smS")
        c_lin["NP", "In-sample", k] <- intAUC(auc_sort_in_lin$NP, 
                                              auc_sort_in_lin$time, 
                                              KM_est_train, method = "smS")
        c_lin["NP", "Out-of-sample", k] <- intAUC_appl(auc_sort_out_lin$NP, 
                                                       auc_sort_out_lin$time, 
                                                       KM_est_test, 
                                                       method = "smS")
        
        ## Smoothed non-parametric
        c_gam["SNP", "In-sample", k] <- intAUC(auc_sort_in_gam$SNP,
                                              auc_sort_in_gam$time, 
                                              KM_est_train, method = "smS")
        c_gam["SNP", "Out-of-sample", k] <- intAUC_appl(auc_sort_out_gam$SNP,
                                                        auc_sort_out_gam$time, 
                                                        KM_est_test, 
                                                        method = "smS")
        c_lin["SNP", "In-sample", k] <- intAUC(auc_sort_in_lin$SNP,
                                               auc_sort_in_lin$time, 
                                               KM_est_train, method = "smS")
        c_lin["SNP", "Out-of-sample", k] <- intAUC_appl(auc_sort_out_lin$SNP,
                                                        auc_sort_out_lin$time, 
                                                        KM_est_test, 
                                                        method = "smS")
        }
       
        setTxtProgressBar(pb, k)
}

##### Figures #####
## concordance
df_c_gam <- as.data.frame.table(c_gam) %>% mutate(model="GAM")
df_c_lin <- as.data.frame.table(c_lin) %>% mutate(model = "Lin")

bind_rows(df_c_gam, df_c_lin) %>% 
  mutate(type = ifelse(estimator %in% c("Harrell", "NP","SNP"), 
                       "Non-parametric", "Semi-parametric")) %>% 
  mutate(model = factor(model, levels=c("GAM", "Lin"), labels = c("ACM", "LCM"))) %>%
        ggplot() + 
        geom_boxplot(aes(x=estimator,y=Freq, fill=type)) + 
        facet_grid(rows = vars(model), cols = vars(estimand)) +
        xlab("") + 
        ylab("Concordance") + labs(fill="")+
  scale_fill_manual(values = cbPalette)+
  theme(text = element_text(size = 15), axis.text = element_text(size = 10))
ggsave(filename = "concordance.png", path = "Images/data_appl", width=10, height=8, bg = "white")

# check concordance numbers
## gam 
df_c_gam %>% group_by(estimator, estimand) %>% summarise_at("Freq", mean)


#bind_rows(df_c_gam, df_c_lin) %>% 
#  group_by(estimator, estimand, model) %>% 
#  summarise_at("Freq", mean)
  

## TV-AUC
tvauc_in_df_gam <- lapply(tvauc_in_gam, as.data.frame) %>%
  bind_rows(.id = "iter") %>%
  mutate(sample = "In-sample", model = "GAM")

tvauc_out_df_gam <- lapply(tvauc_out_gam, as.data.frame) %>%
  bind_rows(.id = "iter") %>%
  mutate(sample = "Out-of-sample" ,model = "GAM")

tvauc_in_df_lin <- lapply(tvauc_in_lin, as.data.frame) %>%
  bind_rows(.id = "iter") %>%
  mutate(sample = "In-sample", model = "Lin")

tvauc_out_df_lin <- lapply(tvauc_out_lin, as.data.frame) %>%
  bind_rows(.id = "iter") %>%
  mutate(sample = "Out-of-sample" ,model = "Lin")

bind_rows(tvauc_in_df_gam, tvauc_out_df_gam, tvauc_in_df_lin, tvauc_out_df_lin) %>% 
  mutate(model = factor(model, levels = c("GAM", "Lin"), labels = c("ACM","LCM"))) %>%
  pivot_longer(c("HZ", "NP", "SNP"), names_to = "Estimator", values_to = "AUC") %>%
  ggplot(aes(x = time, y = AUC, col = model, linetype = sample))+
  geom_smooth(se = F, formula = y~s(x,  k=30, bs = "cs"),
              na.rm = T, method = "gam")+
  labs(x = "Time", y = "AUC")+
  facet_grid(~Estimator)+
  scale_color_manual(values = cbPalette)+
  theme(text = element_text(size = 15), axis.text = element_text(size = 10))
ggsave(filename = "tvauc_gam.png", path = "Images/data_appl", 
       width=15, height=4, bg="white")


save(df_c_gam, df_c_lin, tvauc_in_df_gam, tvauc_in_df_lin, tvauc_out_df_gam, tvauc_out_df_lin,
     file = here("Data/ApplResult.RData"))
# out-of-sample SNP of ACM seems wired (too high)

## only one fold
## fold 1, 5  are abnormally high
# many folds seems to have concordance closer to later (higher values)
tvauc_out_df_gam %>% select(time, SNP, model, iter) %>% filter(model == "GAM") %>%
  ggplot(aes(x = time, y = SNP))+
  geom_point()+
  facet_wrap(~iter)

tvauc_in_df_gam %>% select(time, SNP, model, iter) %>% filter(model == "GAM") %>%
  ggplot(aes(x = time, y = SNP))+
  geom_point()+
  facet_wrap(~iter)

df_c_gam %>% filter(estimator=="SNP" & estimand == "Out-of-sample")

bind_rows(tvauc_out_df_gam %>% select(time, SNP, model, iter) %>% filter(model == "GAM") %>% mutate(sample="out"),
          tvauc_in_df_gam %>% select(time, SNP, model, iter) %>% filter(model == "GAM") %>% mutate(sample="in")) %>%
  ggplot(aes(x = time, y = SNP, col=sample))+
  geom_point()+
  facet_wrap(~iter)
  

## suvival functions
par(mfrow = c(2, 5))
for(k in 1:nfolds){

  df_test_k <- df_analysis_subj[inx_ls[[k]],]
  KM_fit_test  <- survfit(Surv(event_time_years,mortstat) ~ 1, timefix=FALSE, data=df_test_k)
  KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]
  t_uni_test <- unique(df_test_k$event_time_years[df_test_k$mortstat==1])
  fit_S <- scam(St ~ s(ut, bs="mpd", k = 20),
                data=data.frame(ut=t_uni_test, St=KM_est_test))
  df_pred <- data.frame(ut = rep(t_uni_test, each=2) + rep(c(-0.00001,0.00001), length(t_uni_test)))
  # numerator
  St_pred <- predict(fit_S, newdata=df_pred, type='response')
  ft     <- -diff(St_pred)[seq(1,2*length(t_uni_test),by=2)]   
  wt <- 2 * ft * KM_est_test
  # denominator
  
  
  plot(t_uni_test, KM_est_test, main = k)
  
}

