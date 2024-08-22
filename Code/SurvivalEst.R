# This script explores the influnce of basis function
# on the smoothed survival function

#### Set up ####

rm(list = ls())
library("survival")
library("ggplot2")
library(ggpubr)
library("cubature")
library(tidyverse)
library(dplyr)
library(here)
library(risksetROC)
library(mgcv)
library(scam)
library("clinfun")
library(glmnet)
library(kableExtra)
theme_set(theme_minimal())
cbPalette <- c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set the seed
set.seed(512)

load(here("Data/SimNoiseTrainData.RData"))

df <- train_list[[1]] 

#### Survival curves ####

# KP
kp_surv <- survfit(Surv(time,event)~1, timefix=FALSE, data=df)
table(df$event)
plot(kp_surv)
kp_surv$surv[kp_surv$n.event>0]

df_surv <- data.frame(time = kp_surv$time[kp_surv$n.event>0],
                      kp = kp_surv$surv[kp_surv$n.event>0])

df_surv2 <- expand_grid(df_surv, 
                       k = c(10, 20, 30),
                       basis = c("MPD", "MDCX"))

# SCAM smoothing
# shape constrained smoothing terms (SCOP-splines)

k_vec <- c(10, 20, 30)
df_surv2$SmSurv <- NA

for(i in seq_along(k_vec)){
  # MPD
  mpd_surv <- scam(kp ~ s(time, bs="mpd",k=k_vec[i]),
                    data=df_surv)
  df_surv2$SmSurv[df_surv2$basis=="MPD"&df_surv2$k==k_vec[i]] <- predict(mpd_surv)
  
  # MPCX
  mpcx_surv <- scam(kp ~ s(time, bs="mdcx",k=k_vec[i]),
                   data=df_surv)
  df_surv2$SmSurv[df_surv2$basis=="MDCX"&df_surv2$k==k_vec[i]] <- predict(mpcx_surv)
  
  
}



# Figure

df_surv2 %>% 
  mutate(k = paste0("k = ", k)) %>%
  ggplot()+
  geom_step(aes(x=time, y=kp, col = "Kaplan-Meier"))+
  geom_line(aes(x=time, y=SmSurv, col = "Smooth"))+
  labs(y = "", x="Time", col = "")+
  scale_color_manual(values = cbPalette)+
  facet_grid(rows = vars(basis), cols = vars(k))
ggsave(filename = here("Code/SurvSmooth.pdf"), 
       width=12, height=8, bg="white", dpi = 300)
  


