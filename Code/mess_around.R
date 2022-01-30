AUC <- auc_sort_test$HZ_eta1
utimes <- auc_sort_test$time
St <- KM_est_test
plot(utimes, St)
points(auc_sort_train$time, KM_est_train, col = "red")
