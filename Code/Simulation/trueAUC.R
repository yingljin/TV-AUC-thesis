# This script calculates true I/D AUC based on data generating mechanism
# which is Weibull(2, 2)

#### functions for calculating true I/D AUC ####

# Need two different functions since there is a double integral in the specificity formula
# x is linear predictor, t is survival time
# lambda and p are scale and shape parameters of Weibull distribution
# mean_eta is the mean of linear predictor
# sigma_eta is the standard error of linear predictor, generated from N(0,1)

# function for integrating incident sensitivity
my_fn_sens <- function(x, t, lambda=2, p=2, mn_eta=0, sigma_eta){
  lambda*p*t^(p-1)*exp(x[1])*exp(-lambda*t^p*exp(x[1]))*dnorm(x[1], mean=mn_eta, sd=sigma_eta)
}
## function for integrating dynamic specificity
my_fn_spec <- function(x, lambda=2, p=2, mn_eta=0, sigma_eta){
  lambda*p*x[2]^(p-1)*exp(x[1])*exp(-lambda*x[2]^p*exp(x[1]))*dnorm(x[1], mean=mn_eta, sd=sigma_eta)
}


##### functions for calculate true survival #####

### marginal survival function
true_st <- function(t, lambda=2, p=2, eta, sigma_eta){
  St <- exp(-(lambda*exp(eta)*t)^p)*dnorm(eta, 0, sigma_eta)
  return(St)
}

### marginal density function

true_ft <- function(t, lambda=2, p=2, eta, sigma_eta){
  ft <- p*lambda^p*t^(p-1)*exp(p*eta-(lambda*exp(eta)*t)^p)*dnorm(eta, 0, sigma_eta)
  return(ft)
}

#### functions for calculate AUC by integrating the ROC curve ####
trap_integrate_ROC <- function(eta, sens, spec){
  if(!Inf %in% eta){
    sens <- c(sens, 0)
    spec <- c(spec, 1)
  }
  if(!c(-Inf) %in% eta){
    sens <- c(1, sens)
    spec <- c(0, spec)
  }
  x <- rev(1-spec)
  y <- rev(sens)
  n <- length(y)
  sum( (x[2:n] - x[1:(n-1)] ) * (y[1:(n-1)] + y[2:n]) )/2
}


##### Calcualte true AUC #####

tind <- seq(0, 1, length.out = 500)
etaind <- seq(-5, 5, length.out = 100)
true_auc <- rep(NA, length(tind))
Beta <- c(1,-1,0.25)

## in calculating sensitivity and specificity
## restrict boundaries of eta and t to a reasonable value based on distribution
## to avoid numeric issues
pb <- txtProgressBar(min=0, max=length(tind), style=3)
for(i in seq_along(tind)){
  t_i <- tind[i]
  ## containers
  sens <- rep(NA, length(etaind))
  spec <- rep(NA, length(etaind))
  for(j in seq_along(etaind)){
    eta_ij <- etaind[j]
    
    sens[j] <- adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(eta_ij), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral/
      adaptIntegrate(my_fn_sens, t=t_i, lowerLimit=c(-Inf), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral
    
    spec[j] <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(eta_ij, 500), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral/
      adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, t_i), upperLimit=c(500, 500), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral
  }
  ## integrate to AUC
  true_auc[i] <- trap_integrate_ROC(etaind, sens, spec)
  
  setTxtProgressBar(pb, value=i)
}

# brief look at true auc
plot(tind, true_auc) 
true_auc_sort <- data.frame(time_bin = tind, auc = true_auc)

##### Calculate true concordance #####

## marginal survival
true_marg_st <- rep(NA, length(tind))
for(i in seq_along(true_marg_st)){
  true_marg_st[i] <- adaptIntegrate(my_fn_spec, lowerLimit=c(-Inf, tind[i]), upperLimit=c(Inf, 500), lambda=2, p=2, sigma_eta =  sqrt(sum(Beta^2)))$integral
}

plot(tind, true_marg_st)

## marginal density
true_marg_ft <- rep(NA, length(tind))
for(i in seq_along(true_marg_ft)){
  true_marg_ft[i] <- adaptIntegrate(my_fn_sens, t=tind[i], lowerLimit=c(-Inf), upperLimit=c(500), lambda=2, p=2, sigma_eta=sqrt(sum(Beta^2)))$integral
}

plot(tind, true_marg_ft)

## use trapezoidal rule to approximate integral
## use truncated version 
w_vec <- 2 * true_marg_ft * true_marg_st
w_vec <- w_vec/(1-true_marg_st[which.max(tind)]^2)
y_vec <- w_vec*true_auc_sort$auc


nt <- nrow(true_auc_sort)
width <- diff(true_auc_sort$time_bin)
height <- y_vec[1:nt-1]+y_vec[2:nt]
true_c <- sum(width*height/2, na.rm = T)

## save results

save(true_auc_sort, true_c, true_marg_ft, true_marg_ft, file = here("outputData/true_values.RData"))
