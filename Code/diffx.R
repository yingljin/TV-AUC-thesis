

# one iteration

hist(data_test$time, breaks = 30)
hist(data_test2$time, breaks = 30)
hist(data_test3$time, breaks = 30)

boxplot(true_eta~event, data_test3)

auc_lst_test_diffx[[1]]


test_auc_est_truebeta <- tv_auc(data_test3$X %*% Beta, data = data_test3, t=t_uni_test3, nt=nt_uni_test3)
test_auc_est5 <- tv_auc(data_test3$X %*% coef(fit_1), data_test3, t_uni_test3, nt_uni_test3)

ggplot(data.frame(test_auc_est_truebeta))+
  geom_line(aes(x=time, y=HZ, col = "HZ"))+
  geom_line(aes(x=time, y=NP, col = "NP"))+
  geom_line(aes(x=time, y=SNP, col = "SNP"))

test_concord4 <- concord(data_test3, KM_est_test3, test_auc_est_truebeta, fit_1$coefficients, data_test3$X)
test_concord5 <- concord(data_test3, KM_est_test3, test_auc_est5, fit_1$coefficients, data_test3$X)  
  
