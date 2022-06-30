
# likelihood function based on survival data
# standard deviation all set to 1
llh <- function(t = df$time, m = df$M, rho = -0.7, cr = 0.4){
  # calculate mean of censor time 
  mu_censor <- -qnorm(cr)*sqrt(1+1)
  # covaraice matrix of bivariate normal
  sigma_mat <- matrix(c(1, rho, rho, 1), byrow = T, nrow = 2)
  # pdf
  pdf <- rep(NA ,N)
  for(i in 1:N){
    pdf[i] <- (pmvnorm(lower = c(log(t[i]), m[i]), upper = c(Inf, m[i]), sigma = sigma_mat)*dnorm(log(t[i]), mu_censor, 1)+
                 dmvnorm(c(log(t[i]), m[i]),sigma = sigma_mat)*(1-pnorm(log(t[i]), mu_censor, 1)))
  }
  
  return(sum(log(pdf)))
}

# test function
llh(df$time, df$M)

# get maximum of rho and then
# caclulate estimated AUC accordingly
mle_rho <- optimize(f = llh, interval = c(0, 1), t = df$time, m = df$M, cr = 0.4, maximum = T)$maximum
ut <- unique(df$time[df$status==1]) # unique event time
c_vec <-  unique(df$M) # unique biomarker values as cutoff sequence
mle_auc <- rep(NA, length(ut))
for(i in seq_along(ut)){
  logt <- log(ut[i])
  # sensitivity: formula from HZ paper in Section 2.5
  true_tp <- pnorm((mle_rho*logt-c_vec)/sqrt(1-mle_rho^2))
  # 1-specificity: formula from HZ paper in Section 2.5
  mat <- matrix(c(c_vec, rep(logt,length(c_vec))), ncol = 2)
  true_fp <- apply(mat, 1, pmvnorm, upper = Inf, mean = c(0, 0),
                   corr = matrix(c(1, mle_rho, mle_rho, 1), nrow = 2, ncol = 2))/pnorm(-logt)
  # AUC
  id <- order(c_vec)
  true_tp <- true_tp[id]
  true_fp <- true_fp[id]
  dFP <- abs(true_fp[-1] - true_fp[-length(true_fp)])
  aTP <- 0.5 * (true_tp[-1] + true_tp[-length(true_tp)])
  mle_auc[i] <- sum(dFP * aTP)
}

