#### plot residual ####
auc_lst[[1]]

par(mfrow = c(3, 3))
mtext("Residual vs Fitted for regular GAM", outer = TRUE)
for(i in 1:9){
  df <- as.data.frame(auc_lst[[i]])
  sm_fit <- gam(empirical ~ s(time, k = 30), data = df)
  pred <- predict(sm_fit, type = "response", se.fit = T)
  plot(pred$fit, sm_fit$residuals, cex = 0.5, pch = 20)
  abline(h = 0)
}

i <- 1

##### random code #####

auc_lst[2:10] %>%
  bind_rows(.id = "iter") %>%
  mutate(true = approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, xout = time)$y) %>%
  pivot_longer(3:6, names_to = "estimator", values_to = "AUC") %>%
  filter(estimator != "empirical") %>%
  ggplot(aes(x = time, y = AUC, group = estimator, col = estimator))+
  geom_line()+
  facet_wrap(~iter)


for(i in 2:10){
  df <- auc_lst[[i]]
  df$empirical<-(1-df$empirical)^(1/3)
  sm_fit1 <- gam(empirical ~ s(time, k = 30), data = df)
  pred_y <- predict(sm_fit1, type = "response", se.fit = F)
  pred_y <- 1-pred_y^3
  auc_lst[[i]]$sm_empirical_trans <- pred_y
}

plot(predict(sm_fit1), sm_fit1$residuals)
hist(sm_fit1$residuals)
hist((1-auc_lst[[1]]$empirical)^(1/3))
hist(auc_lst[[100]]$empirical)

head(auc_lst[[i]])

auc_lst[2:10] %>%
  bind_rows(.id = "iter") %>%
  mutate(true = approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, xout = time)$y)%>%
  pivot_longer(3:8, names_to = "estimator", values_to = "AUC") %>%
  filter(estimator %in% c("sm_empirical", "sm_empirical_trans", "true")) %>%
  ggplot(aes(x = time, y = AUC, group = estimator, col = estimator))+
  geom_line()+
  facet_wrap(~iter)

auc_lst[2:10] %>%
  bind_rows(.id = "iter") %>%
  mutate(true = approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, xout = time)$y) %>%
  mutate(bias_hz = HZ-true,
         bias_sm_em = sm_empirical-true,
         bias_sm_em_weight = sm_empirical_weighted-true) %>%
  dplyr::select(iter, time, bias_hz, bias_sm_em, bias_sm_em_weight) %>%
  pivot_longer(3:5, names_to = "estimator", values_to = "bias") %>%
  ggplot(aes(x=time, y=bias, group = estimator, col = estimator))+
  geom_point()
