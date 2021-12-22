
library(tidyverse)
library(mgcv)
library(here)

load("outputData/results_by_iter.RData")
load("outputData/true_values.RData")

#### plot residual ####
auc_lst[[1]]

par(mfrow = c(3, 3))
for(i in 1:9){
  df <- as.data.frame(auc_lst[[i]])
  sm_fit <- gam(empirical ~ s(time, k = 30), data = df)
  pred <- predict(sm_fit, type = "response", se.fit = T)
  plot(pred$fit, sm_fit$residuals, cex = 0.5, pch = 20)
  abline(h = 0)
}
mtext("Residual vs Fitted for unweighted GAM", outer = TRUE, cex = 0.8, line = -2)

##### weighted gam #####

# calculate variance for weights
# use empirical estimators
## bin into small bins 
auc_df <- auc_lst %>% bind_rows(.id = "iter")
head(auc_df)
auc_df <- auc_df %>% 
  mutate(time_bin = cut(time, breaks=seq(0,1,length.out = 101), labels = 1:100, include.lowest = T))
table(auc_df$time_bin) 
## calculate variance in each bin
var_df <- auc_df %>%
  group_by(time_bin) %>%
  summarise(bin_var = var(empirical)) %>%
  mutate(wt = 1/bin_var)
hist(var_df$bin_var, breaks = 15)
hist(var_df$wt, breaks = 15)
## put time and corresponding weight together
wt_df <- auc_df %>% 
  select(time, time_bin) %>%
  left_join(var_df, by = "time_bin")

# Weighted gam
for(i in seq_along(auc_lst)){
  df <- as.data.frame(auc_lst[[i]]) %>%
    left_join(wt_df, by = "time")
  sm_fit <- gam(empirical ~ s(time, k = 30), data = df, weights = wt)
  pred <- predict(sm_fit, type = "response", se.fit = T)
  auc_lst[[i]] <- 
    auc_lst[[i]] %>% mutate(sm_empirical_wt = pred$fit)
}

# export estimators
save(auc_lst, file = here("outputData/weighted_gam_by_iter.RData"))
