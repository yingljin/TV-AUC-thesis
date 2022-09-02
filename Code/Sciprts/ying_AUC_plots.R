rm(list=ls())
library(ggplot2)
library(tidyverse)
load("~/Downloads/TVAUCestimates.RData")


plt_in <- 
        tvauc_in_df %>% 
        pivot_longer(cols=c("HZ","NP","SNP")) %>% 
        ggplot() + 
        geom_point(aes(x=time, y=value)) + 
        facet_grid(~name) + 
        geom_smooth(aes(x=time, y=value), method="gam", formula=y~s(x,bs="cr",k=20))

plt_out <- 
        tvauc_out_df %>% 
        pivot_longer(cols=c("HZ","NP","SNP")) %>% 
        ggplot() + 
        geom_point(aes(x=time, y=value)) + 
        facet_grid(~name) + 
        geom_smooth(aes(x=time, y=value), method="gam", formula=y~s(x,bs="cr",k=20))


library(mgcv)

fit_np_GCV <- gam(NP ~ s(time, bs="cr",k=20), method="GCV.Cp", data=tvauc_in_df)
fit_np_REML <- gam(NP ~ s(time, bs="cr",k=20), method="REML", data=tvauc_in_df)


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
