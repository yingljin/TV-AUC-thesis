
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

par(mfrow = c(3, 3))
for(i in 1:9){
  df <- as.data.frame(auc_lst[[i]])
  sm_fit <- gam(empirical ~ s(time, k = 30), data = df)
  pred <- predict(sm_fit, type = "response", se.fit = T)
  weights <- 1/pred$se.fit^2
  sm_fit_weighted <- gam(empirical ~ s(time, k = 30), data = df, weights = weights)
  pred_weighted <- predict(sm_fit_weighted, type = "response", se.fit = T)
  plot(pred_weighted$fit, sm_fit_weighted$residuals, cex = 0.5, pch = 20)
  abline(h = 0)
  auc_lst[[i]] <- auc_lst[[i]] %>% mutate(weighted_sm_empirical = pred_weighted$fit)
}
mtext("Residual vs Fitted for weighted GAM", outer = TRUE, cex = 0.8, line = -2)

auc_lst[1:9] %>%
  bind_rows(.id = "iter") %>%
  pivot_longer(3:6, names_to = "estimator", values_to = "AUC") %>%
  filter(estimator %in% c("sm_empirical", "weighted_sm_empirical", "HZ")) %>%
  ggplot(aes(x = time, y = AUC, group = estimator, col = estimator))+
  geom_line()+
  facet_wrap(~iter)+
  labs(title = "Compare smoothed empirical estimator before and after weighting")

##### transformed gam #####

par(mfrow = c(3, 3))
for(i in 1:9){
  df <- as.data.frame(auc_lst[[i]])
  df$trans_em <- log(1-df$empirical)
  df$trans_em[df$trans_em==-Inf] <- NA
  mis_id <- which(!is.na(df$trans_em))
  sm_fit_trans <- gam(trans_em ~ s(time, k = 30), data = df)
  pred_trans <- predict(sm_fit_trans, type = "response", se.fit = T)
  plot(pred_trans$fit[mis_id], sm_fit_trans$residuals, cex = 0.5, pch = 20)
  abline(h = 0)
  auc_lst[[i]] <- auc_lst[[i]] %>% mutate(trans_sm_empirical = NA)
  auc_lst[[i]][mis_id, "trans_sm_empirical"] <- 1-exp(pred_trans$fit)
}
mtext("Residual vs Fitted for transformed GAM", outer = TRUE, cex = 0.8, line = -2)

auc_lst[1:9] %>%
  bind_rows(.id = "iter") %>%
  pivot_longer(3:7, names_to = "estimator", values_to = "AUC") %>%
  filter(estimator %in% c("sm_empirical", "trans_sm_empirical", "HZ")) %>%
  ggplot(aes(x = time, y = AUC, group = estimator, col = estimator))+
  geom_line()+
  facet_wrap(~iter)+
  labs(title = "Compare smoothed empirical estimator before and after transformation")

##### Other random codes #####
par(mfrow = c(1, 1))
hist(true_auc_sort$auc, breaks = 20)

hist(auc_lst[[1]]$empirical, breaks = 100)
