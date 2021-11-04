library(MASS)
library(risksetROC)
library(survival)

# investigaete abnormal value of true auc when when t = 0.145 and eta = 1.76767677
t_i <- 0.145
# t_i <- 0.155
etaind <- seq(-7, 7, length.out = 100)
sens <- rep(NA, length(etaind))
spec <- rep(NA, length(etaind))
for(j in seq_along(etaind)){
  eta_ij <- etaind[j]
  
  sens[j] <- adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(eta_ij), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral/
    adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(-Inf), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral
  
  spec[j] <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(eta_ij, Inf), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral/
    adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(Inf, Inf), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral
}

plot(etaind, sens)
plot(etaind, spec)

# weird specificity when eta = 1.77
eta_ij <- 1.76767677
# eta_ij <- 1.62626263
# eta_ij <- 1.90909091

## specificity
num2 <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(eta_ij, Inf), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))
dem2 <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(Inf, Inf), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))

num2$integral/dem2$integral # 1.013>1

## function to calculate specificity
## x is actually a vector(eta, t)
my_fn_spec <- function(x, lambda=1, p=1, mn_eta=0, sigma_eta=2){
  lambda*p*x[2]^(p-1)*exp(x[1])*exp(-lambda*x[2]^p*exp(x[1]))*dnorm(x[1], mean=mn_eta, sd=sigma_eta)
}

## fix t
eta_vec <- seq(min(etaind), eta_ij, length.out = 1000)
df <- cbind(eta_vec, t_i)
spec_vec <- apply(df, 1, my_fn_spec, lambda = 2, p = 2, sigma_eta = sqrt(sum(Beta^2)))
plot(eta_vec, spec_vec)

## fix eta
t_vec <- seq(t_i, max(tind), length.out = 1000)
df <- cbind(eta_ij, t_vec)
spec_vec <- apply(df, 1, my_fn_spec, lambda = 2, p = 2, sigma_eta = sqrt(sum(Beta^2)))
plot(t_vec, spec_vec)
