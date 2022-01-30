# this script generates motivational example for presentation

##### Linear regression #####

x1 <- rnorm(100)
y1 <- 0.25*x1
y2 <- 0.25*x1+rnorm(100, mean= 0, sd = 0.25)

plot(x1, y1, xlab = "x", ylab = "y")
plot(x1, y2, xlab = "x", ylab = "y")
abline(lm(y2~x1), col = "red")

##### logistic regression #####

x1 <- c(rnorm(50, mean = 0.25), rnorm(50, mean = 0.75))
y1 <- rep(c(0, 1), each = 50) 
lg_fit <- glm(y1~x1, family = binomial)

boxplot(x1~y1, xlab = "Y", ylab = "X")

##### Survival analysis #####

library(survival)

df <- data.frame(id = 1:7, 
                 start = runif(7))

df <- df %>% 
  mutate(end = runif(7, start, 1), 
         status = sample(c("Censor", "Event"), size = 7, prob = c(0.5, 0.5), replace = T))
df$id <- factor(df$id)
ggplot(df)+
  geom_linerange(aes(xmin=start, xmax=end, y=id, colour = status), size = 1.2, show.legend = F)+
  geom_point(aes(x=end, y=id, shape = status, colour = status), size = 5)+
  labs(x = "Time", y = "Observations", title = "Time-to-Event data example")

#### TV-AUC distribution ####
head(auc_df)

auc_df %>% select(time, contains("HZ_")) %>%
  pivot_longer(2:4) %>%
  ggplot(aes(x = value))+
  geom_histogram()+
  facet_wrap(~name)
