
rm(list = ls())
library("survival")
library("ggplot2")
library(ggpubr)
library("cubature")
library("tidyverse")
library(here)
library(risksetROC)
library(mgcv)
library(scam)
library("clinfun")
library(glmnet)
theme_set(theme_minimal())
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set the seed
set.seed(726)



load(here("Data/results_by_iter.RData"))

##### Cleaning and formatting results #####

## AUC
auc_df_train <- bind_rows(auc_lst_train, .id = "iter")
auc_lst_train[[1]] %>% filter(model=="No noise") %>%
  ggplot(aes(x=time, y=HZ))+
  geom_line()

auc_df_test <- bind_rows(auc_lst_test, .id = "iter")
auc_df <- bind_rows(auc_df_train, auc_df_test, .id = "sample")


table(auc_df$model)

## concordance
c_df_train <- bind_rows(c_lst_train, .id = "iter")
c_df_test <- bind_rows(c_lst_test, .id = "iter")
c_df <- bind_rows(c_df_train, c_df_test, .id = "sample")

table(c_df$model)

#### true AUC and concordance #####
load(here("Data/true_values.RData"))

# use interpolation to estimate AUC on simulated time series
true_auc_sort <- approx(x = true_auc_sort$time_bin, y = true_auc_sort$auc, 
                        xout = auc_df$time)$y

head(true_auc_sort)

#### Produce figures #####

# trend of TV-AUC estimates (smoothed)
auc_df_long <- auc_df %>% 
  mutate(true = true_auc_sort,
         sample = factor(sample, levels = 1:2, labels = c("In-sample", "Out-of-sample"))) %>% 
  pivot_longer(4:6) 

head(auc_df_long)

auc_df_long$model <- factor(auc_df_long$model, levels =  c("No noise", "20 noise", "20 noise, penalized", 
                                                           "100 noise", "100 noise, penalized", "Different covariate"))
table(auc_df_long$model)
auc_df_long$name <- factor(auc_df_long$name, levels = c("HZ", "SNP", "NP"))
table(auc_df_long$name)

#### TV-AUC all iterations ####

auc_df_long %>% 
  filter(model!="20 noise, penalized" & model!="100 noise, penalized") %>%
  ggplot(aes(x=time, y=value, col=model, linetype = sample))+
  geom_smooth(se = F, formula = y~s(x, k=30, bs = "cs"), na.rm = T,
              method = "gam")+
  geom_line(aes(x = time, y = true), na.rm = T, col = "red", show.legend = F)+
  labs(x="time", y = "AUC")+
  theme(text = element_text(size = 15), axis.text = element_text(size = 10))+
  facet_wrap(~name)+
  scale_colour_manual(values=cbPalette)

ggsave(filename = "tvauc_N250.png", path = "Images/N250", width=15, height=4, bg = "white")
#ggsave(filename = "tvauc_N500.png", path = "Images/N500", width=15, height=4, bg = "white")

# spread of TV-AUC estimates
brk <- seq(0, 1, 0.2) # bin time into five different bins
ggarrange(
  auc_df_long %>%
    # dplyr::select(-true) %>%
    filter(sample == "In-sample") %>%
    filter(model!="20 noise, penalized" & model!="100 noise, penalized") %>%
    mutate(time_bin = cut(time, breaks = brk, include.lowest = T)) %>%
    ggplot(aes(x = factor(time_bin), y = value, fill = name))+
    geom_boxplot(outlier.size = 0.5)+
    facet_grid(cols = vars(model))+
    labs(x = "Time", y = "AUC", title = "In-sample")+
    theme(axis.text.x = element_text(angle = 60, vjust = 0.1, hjust = 0.1), 
          text=element_text(size = 15),
          axis.text = element_text(size = 10))+
    scale_fill_manual(values=cbPalette)+
    labs(y="AUC"), 
  
  auc_df_long %>%
    # dplyr::select(-true) %>%
    filter(sample == "Out-of-sample") %>%
    filter(model!="20 noise, penalized" & model!="100 noise, penalized") %>%
    mutate(time_bin = cut(time, breaks = brk, include.lowest = T)) %>%
    ggplot(aes(x = factor(time_bin), y = value, fill = name))+
    geom_boxplot(outlier.size = 0.5)+
    facet_grid(cols = vars(model))+
    labs(x = "Time", y = "AUC", title = "Out-of-sample")+
    theme(axis.text.x = element_text(angle = 60, vjust = 0.1, hjust = 0.1), 
          text = element_text(size=15),
          axis.text = element_text(size = 10))+
    scale_fill_manual(values=cbPalette)+
    labs(y="AUC"), nrow=1, common.legend = T, widths = c(3, 4))
ggsave(filename = "tvauc_box_N250.png", path = "Images/N250", width=15, height = 4, bg = "white")
#ggsave(filename = "tvauc_box_N500.png", path = "Images/N500", width=15, height = 4, bg = "white")

# concordance
c_df <- bind_rows(c_df_train, c_df_test, .id = "sample")
c_df_long <- c_df %>%
  mutate(sample = factor(sample, levels = 1:2, labels = c("In-sample", "Out-of-sample"))) %>% 
  pivot_longer(3:7) %>%
  mutate(Type = ifelse(name=="HZ"|name=="GH", "Semi-parametrc", "Non-parametric")) %>%
  mutate(model = factor(model, levels =  c("No noise", "20 noise", "20 noise, penalized", "100 noise", "100 noise, penalized", "Different covariate")),
         name = factor(name, levels = c("HZ","GH", "Harrell.s","NP","SNP"),
                       labels = c("HZ","GH", "Harrell","NP","SNP")))

ggarrange(
  c_df_long %>% 
    filter(sample == "In-sample") %>%
    filter(model!="20 noise, penalized" & model!="100 noise, penalized") %>%
    ggplot(aes(x = name, y = value))+
    geom_boxplot(aes(fill = Type))+
    facet_grid(cols=vars(model))+
    geom_hline(yintercept = true_c, col = "red")+
    theme(text = element_text(size = 12),
          axis.text = element_text(size=8))+
    scale_fill_manual(values=cbPalette)+
    labs(y = "Concordance", x = "Estimator", title = "In-sample"),
  
  c_df_long %>% 
    filter(sample == "Out-of-sample") %>%
    filter(model!="20 noise, penalized" & model!="100 noise, penalized") %>%
    ggplot(aes(x = name, y = value))+
    geom_boxplot(aes(fill = Type))+
    facet_grid(cols=vars(model))+
    geom_hline(yintercept = true_c, col = "red")+
    theme(text = element_text(size = 12),
          axis.text = element_text(size=8))+
    scale_fill_manual(values=cbPalette)+
    labs(y = "Concordance", x = "Estimator", title = "Out-of-sample"), nrow = 1, common.legend = T, widths = c(3, 4))
ggsave(filename = "concordance_N250.png", path = "Images/N250", width=15, height = 4, bg = "white")
#ggsave(filename = "concordance_N500.png", path = "Images/N500", width=15, height = 4, bg = "white")

#### TV-AUC example ####
library(MASS)
data(VA)
survival.time <- VA$stime
survival.status <- VA$status
score <- VA$Karn
cell.type <- factor(VA$cell )
tx <- as.integer( VA$treat==1 )
age <- VA$age
survival.status[VA$stime > 500 ] <- 0
survival.time[VA$stime > 500 ] <- 500
library(survival)
fit0 <- coxph( Surv(survival.time,survival.status)
               ~ score + cell.type + tx + age, na.action=na.omit )
summary(fit0)
eta <- fit0$linear.predictor
AUC <- NULL

event_t <- survival.time[survival.status==1]
event_t <- sort(event_t)
auc_t <- rep(NA, length(event_t))

for(i in seq_along(event_t)){
out <- CoxWeights(marker=eta, Stime=survival.time, status=survival.status,
                  predict.time=event_t[i])
## to see how well the marker predicts one-month survival
auc_t[i] <- out$AUC
}



plot(event_t, auc_t, type = "b", xlab = "Time", ylab = "AUC")


