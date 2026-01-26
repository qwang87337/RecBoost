# default values:
# c = Cvec
# Cvec_lt = 0
# delta = atrisk$nev_id,
# alpha = pars$alpha,
# bbeta = pars$bbeta,
# pvfm = pvfm,
# dist = pars$dist


# Notes:
# 1. for gamma, alpha = beta = theta, the exponential of the estimated theta from outer problem
#
# 2. delta = atrisk$nev_id, which is the number of event per ID(cluster).
#
# 3. the first column and the second column are the numerators and the denominators of the frailty 
# fraction, the third is the log-likelihood contribution, and the last column is the expectation 
# of the squared frailty(only used in calculating the information matrix)
#
# 4. dist_to_pars() looks like this:
# pars <- dist_to_pars(dist, logfrailtypar, pvfm)
# dist_to_pars <- function(dist, logfrailtypar, pvfm) {
#  
#  if (dist == "gamma") {
#    alpha <- bbeta <- exp(logfrailtypar)
#    dist_id <- 0L
#  }
#}

fast_Estep <- function(c, c_lt = 0, delta, alpha, bbeta, pvfm, dist) {
  
  
  res <- matrix(0, length(delta), 4)
  
  if(dist==0) {
    bbeta <- bbeta + c_lt
    
    
    # res[,3] is the log of kth derivative of gamma's laplace transform. See 
    # appendix A.2 in page 26 of Balan et al., 2019
    res[,3] <- alpha * log(bbeta) - (alpha + delta)*log(bbeta + c) + lgamma(alpha + delta) - lgamma(alpha)
    res[,1] <- (alpha + delta)
    res[,2] <- (bbeta + c)
    res[,4] <- (alpha + delta) * (alpha + delta + 1) / (bbeta + c)^2
  }
  
  res
}
