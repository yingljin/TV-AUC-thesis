
# See what happen to the abonarmally high out-of-sample SNP TV-AUC from GAM model
# In data application

#### clean results #####
df_c_gam <- as.data.frame.table(c_gam) %>% mutate(model="GAM")
df_c_lin <- as.data.frame.table(c_lin) %>% mutate(model = "Lin")



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

##### Plot by iteration ####
### out-of-sample ACM SNP
df_c_gam %>% filter(estimator=="SNP" & estimand == "Out-of-sample")

tvauc_out_df_gam %>% filter(sample == "Out-of-sample") %>% 
  left_join(df_c_gam %>% filter(estimator=="SNP" & estimand == "Out-of-sample") %>%
               rename("iter" = "fold") %>% select(iter, Freq), by = "iter") %>%
  ggplot()+
  geom_point(aes(x = time, y = SNP))+
  geom_hline(aes(yintercept = Freq, col="red"))+facet_wrap(~iter)


### out-of-sample LCM SNP
df_c_gam %>% filter(estimator=="SNP" & estimand == "In-sample")

tvauc_out_df_lin %>% filter(sample == "Out-of-sample") %>% 
  left_join(df_c_lin %>% filter(estimator=="SNP" & estimand == "Out-of-sample") %>%
              rename("iter" = "fold") %>% select(iter, Freq), by = "iter") %>%
  ggplot()+
  geom_point(aes(x = time, y = SNP))+
  geom_hline(aes(yintercept = Freq, col="red"))+facet_wrap(~iter)

### out-of-sample AUC and weights

# for(k in 1:nfolds){
k <- 1 
df_train_k <- df_analysis_subj[-inx_ls[[k]],]
df_test_k <- df_analysis_subj[inx_ls[[k]],]

t_uni_test <- unique(df_test_k$event_time_years[df_test_k$mortstat==1])
nt_uni_test <- length(t_uni_test)

# KP survival
KM_fit_test  <- survfit(Surv(event_time_years,mortstat) ~ 1, timefix=FALSE, data=df_test_k)
KM_est_test <- KM_fit_test$surv[KM_fit_test$n.event>0]

plot(KM_fit_test)
summary(KM_fit_test, time = 1)$surv

# function to calculated weight
tao <- summary(KM_fit_test, time = 1)$surv
weighted_auc <- function(AUC, utimes, St){
  ut <- utimes
  
  # smooth survival function
  fit_S <- scam(St ~ s(ut, bs="mpd", k = 20),
                data=data.frame(ut=ut, St=St))
  df_pred <- data.frame(ut = rep(ut, each=2) + rep(c(-0.00001,0.00001), length(ut)))
  St_pred <- predict(fit_S, newdata=df_pred, type='response')
  St_sm <- predict(fit_S, newdata=data.frame(ut=ut), type='response')
  ft     <- -diff(St_pred)[seq(1,2*length(ut),by=2)]                    

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
  
  return(list(iAUC = iAUC, wt=wt, W=W))
}

tvauc_out_df_gam %>% filter(iter==1)
auc_sort_out_gam <-arrange(tvauc_out_df_gam %>% filter(iter==1), time)
fold1 <- weighted_auc(auc_sort_out_gam$SNP,
            auc_sort_out_gam$time, 
            KM_est_test)



# weight with the KM survival

par(mfrow=c(2, 2))
plot(auc_sort_out_gam$time, auc_sort_out_gam$SNP, main = "out-of-sample SNP", ylab = "")
plot(auc_sort_out_gam$time, fold1$wt*auc_sort_out_gam$SNP/fold1$W, 
     main = "out-of-sample weighted SNP", ylab="")

sum(auc_sort_out_gam$SNP*fold1$wt)


# compare to in-sample
df_train_k <- df_analysis_subj[-inx_ls[[k]],]

t_uni_test <- unique(df_test_k$event_time_years[df_test_k$mortstat==1])
nt_uni_test <- length(t_uni_test)

KM_fit_train <- survfit(Surv(event_time_years,mortstat) ~ 1, timefix=FALSE, data=df_train_k)
KM_est_train <- KM_fit_train$surv[KM_fit_train$n.event>0]

auc_sort_in_gam <-arrange(tvauc_in_df_gam %>% filter(iter==k), time)
fold1_train <- weighted_auc(auc_sort_in_gam$SNP,
                      auc_sort_in_gam$time, 
                      KM_est_train)

plot(auc_sort_in_gam$time, auc_sort_in_gam$SNP, main = "in-sample SNP", ylab="")
plot(auc_sort_in_gam$time, fold1_train$wt*auc_sort_in_gam$SNP/fold1_train$W,
     main = "in-sample weighted SNP", ylab="")
sum(fold1_train$wt*auc_sort_in_gam$SNP/fold1_train$W, main = "out-of=sample SNP")
