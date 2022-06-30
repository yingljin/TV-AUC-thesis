# this script generates motivational example for presentation

##### Linear regression #####

x1 <- rnorm(100)
y1 <- 0.25*x1
y2 <- 0.25*x1+rnorm(100, mean= 0, sd = 0.25)

plot(x1, y1, xlab = "Predictor", ylab = "Outcome")
abline(lm(y1~x1), col = "red")
plot(x1, y2, xlab = "Predictor", ylab = "Outcome")
abline(lm(y2~x1), col = "red")

##### logistic regression #####

x1 <- c(rnorm(50, mean = 0.25), rnorm(50, mean = 0.75))
x2 <- c(rnorm(50, mean = 0.1, sd = 0.1), rnorm(50, mean = 0.9, sd= 0.1))
y1 <- rep(c("Survive", "Death"), each = 50) 
lg_fit <- glm(y1~x1, family = binomial)

boxplot(x1~y1, ylab = "Predictor", xlab = "")
boxplot(x2~y1, ylab = "Predictor", xlab = "")

##### Survival analysis #####

library(survival)
library(tidyverse)

df <- data.frame(id = 1:7,
                 start = rep(0, 7))

df <- df %>% 
  mutate(end = runif(7, start,1), 
         status = sample(c("Censor", "Event"), size = 7, prob = c(0.5, 0.5), replace = T))
df$id <- factor(df$id)
ggplot(df)+
  geom_linerange(aes(xmin=start, xmax=end, y=id, colour = status), size = 1.2, show.legend = F)+
  geom_point(aes(x=end, y=id, shape = status, colour = status), size = 5)+
  labs(x = "Time", y = "Observations", title = "Time-to-Event data example")

#### TV-AUC distribution ####
head(auc_df)

auc_df %>% select(time, true, contains("_eta1")) %>% 
  pivot_longer(2:5) %>%
  mutate(name = factor(name, levels = c("true", "HZ_eta1", "empirical_eta1", "sm_empirical_eta1"),
                       labels = c("True", "HZ", "NP", "SNP"))) %>%
  ggplot(aes(x = value))+
  geom_histogram()+
  facet_wrap(~name)+
  labs(x = "AUC", y = "")


##### some more figures for presentation #####

# TV-AUC
p1
p2
p3

# concordance
ggpubr::ggarrange(pc1, pc2, common.legend = T)
