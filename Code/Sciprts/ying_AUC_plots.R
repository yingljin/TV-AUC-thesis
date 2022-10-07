rm(list=ls())
library(ggplot2)
library(tidyverse)
library(here)
library(arsenal)
list.files("OutputData")
load("OutputData/appl_tvauc.RData")
load("OutputData/appl_c.RData")


#### Concordance #####
df_c_lin %>% 
 group_by(estimand, estimator) %>% 
  summarise_at("Freq", mean)


#### TV-AUC ####
# why is in-sample curve much wigglier than out-of-sample
# The wiggliness is not because of dips? (still wiggly after removing the dips)
# but it is also not because of sample size (still wiggly after random down-sampling)
plt_in <- 
        tvauc_in_df_gam %>% 
        pivot_longer(cols=c("HZ","NP","SNP")) %>% 
        ggplot() + 
        geom_point(aes(x=time, y=value, group = iter, col = iter)) + 
        facet_grid(~name) + 
        geom_smooth(aes(x=time, y=value, group = iter),  se = F, method="gam", formula=y~s(x,bs="cr",k=30))

plt_out <- 
        tvauc_out_df_gam %>% 
  filter(NP > 0.6) %>%
        pivot_longer(cols=c("HZ","NP","SNP")) %>% 
        ggplot() + 
        geom_point(aes(x=time, y=value, group = iter, col = iter)) + 
        facet_grid(~name) + 
        geom_smooth(aes(x=time, y=value),se = F,  method="gam", formula=y~s(x,bs="cr",k=20))


in_NP <- tvauc_in_df_gam %>% select(iter, time, NP)
out_NP <- tvauc_out_df_gam %>% select(iter, time, NP)
bind_rows(in_NP, out_NP, .id = "sample") %>% 
  ggplot()+
  geom_point(aes(x=time, y=NP, col = sample), size = 0.5)

bind_rows(in_NP, out_NP, .id = "sample") %>% 
  ggplot()+
  geom_boxplot(aes(x=time, y=NP, group = time))+
  #geom_jitter(aes(x=time, y=NP))+
  facet_wrap(~sample)

#### Notes ####

# out-of-sample have large variation at all time-points, spread out across iterations
# in-sample estimates have smaller variation (centered across iterations)
# ggsmooth is just captureing a essential trend
# why did that happen? 
# some iterations has much larger variation in outcome



in_HZ <- tvauc_in_df_gam %>% select(iter, time, HZ)
out_HZ <- tvauc_out_df_gam %>% select(iter, time, HZ)
bind_rows(in_HZ, out_HZ, .id = "sample") %>% 
  ggplot()+
  geom_point(aes(x=time, y=HZ, col = sample), size = 0.5)

bind_rows(in_HZ, out_HZ, .id = "sample") %>% 
  ggplot()+
  geom_boxplot(aes(x=time, y=HZ, group = time))+
  #geom_jitter(aes(x=time, y=NP))+
  facet_wrap(~sample)


in_SNP <- tvauc_in_df_gam %>% select(iter, time, SNP)
out_SNP <- tvauc_out_df_gam %>% select(iter, time, SNP)
bind_rows(in_SNP, out_SNP, .id = "sample") %>% 
  ggplot()+
  geom_point(aes(x=time, y=SNP, col = sample), size = 0.5)

bind_rows(in_SNP, out_SNP, .id = "sample") %>% 
  ggplot()+
  geom_boxplot(aes(x=time, y=SNP, group = time))+
  #geom_jitter(aes(x=time, y=NP))+
  facet_wrap(~sample)


#### covariates by time ####
plot(df_analysis_subj$event_time_years, df_analysis_subj$ASTP_mean)
plot(df_analysis_subj$event_time_years, df_analysis_subj$RA_MIMS_mean)
plot(df_analysis_subj$event_time_years, df_analysis_subj$age_years_interview)
plot(df_analysis_subj$event_time_years, df_analysis_subj$BMI)
plot(df_analysis_subj$event_time_years, df_analysis_subj$TMIMS_mean)


#### Other #### 


library(mgcv)

# fit_np_GCV <- gam(NP ~ s(time, bs="cr",k=20), method="GCV.Cp", data=tvauc_in_df)
fit_np_REML <- gam(NP ~ s(time, bs="cr",k=30), method="REML", data=tvauc_in_df_gam)
data.frame(time = tvauc_in_df_gam$time, y = predict(fit_np_REML)) %>% 
  ggplot()+
  geom_boxplot(aes(x=time, y=y, group = time))


k <- 10
fit_ls <- vector(mode="list",length=10)
for(k in 1:10){
        df_k <- subset(tvauc_in_df, iter == k)
        fit_ls[[k]] <- gam(NP ~ s(time, bs="cr",k=20), method="REML", data=df_k)
}


tind_pred <- seq(0,8,len=100)
auc_hat_mat <- matrix(NA, 10, 100)
par(mfrow=c(3,4))
for(k in 1:10){
        plot(fit_ls[[k]], ylim=c(0.5,1), shift=coef(fit_ls[[k]])[1])
        auc_hat_mat[k,] <- predict(fit_ls[[k]],type='response', newdata=data.frame(time=tind_pred))
}

plot(colMeans(auc_hat_mat))


par(mfrow=c(1,2))
plot(fit_np_GCV)
plot(fit_np_REML)
