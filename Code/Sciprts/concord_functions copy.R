## Write a function to calculate Concordance using concordant vs discordant pairs empricical estimator
## Could use survival::concordance, but for some reason I'm getting an error when trying to predict out-of-sample.
calc_c <- function(marker, Stime, status){
  utimes <- sort(unique(Stime[status==1]))
  num <- denom <- 0
  for(ut in seq_along(utimes)){
    ## current time
    ti    <- utimes[ut]
    ## subjects who experienced an event at current time
    inx_i <- which(Stime == ti & status==1)
    ## subjects with observed times beyond event current time
    inx_j <- which(Stime > ti)
    ## number of "cases" and "controls" at time current time
    n_case_t    <- length(inx_i)
    n_control_t <- length(inx_j)
    for(i in seq_along(inx_i)){
      num   <- num + sum( (marker[inx_j] > marker[inx_i[i]] ) ) + 0.5*sum(marker[inx_j] == marker[inx_i[i]])
    }
    denom <- denom + n_case_t*n_control_t
  }
  1-num/denom
}

## calculate concordance under Gonen and heller
calc_c_gh <- function(beta_hat, X){
  N <- nrow(X)
  C <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      eta_ij <- (X[i,,drop=F] - X[j,,drop=F]) %*% beta_hat
      eta_ji <- (X[j,,drop=F] - X[i,,drop=F]) %*% beta_hat
      if(eta_ij < 0){
        C <- C + 1/(1 + exp(eta_ij))
      } 
      if(eta_ji < 0){
        C <- C + 1/(1+exp(eta_ji))
      }
    }
  }
  2/(N*(N-1))*C
}