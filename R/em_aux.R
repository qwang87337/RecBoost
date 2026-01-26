dist_to_pars <- function(dist, logfrailtypar, pvfm) {
  
  if (dist == "gamma") {
    alpha <- bbeta <- exp(logfrailtypar)
    dist_id <- 0L
  }
  
  # if (dist == "stable") {
  #     theta <- exp(logfrailtypar) + 1  # so theta >1
  #     bbeta <- 1 - 1/theta  # so bbeta in (0,1), that's what's important
  #     alpha <- theta / (theta - 1)  # alpha = 1/beta for scaling
  #     dist_id <- 1L
  # }
  
  if (dist == "stable") {
    # theta <- exp(logfrailtypar) + 1 # so theta >1
    # bbeta <- 1 - 1/theta
    alpha <- 1
    #bbeta <- 1 - exp(logfrailtypar) / (exp(logfrailtypar) + 1)
    bbeta <- exp(logfrailtypar) / (exp(logfrailtypar) + 1)
    dist_id <- 1L
  }
  
  if (dist == "pvf") {
    alpha <- abs((pvfm + 1)/pvfm * exp(logfrailtypar))
    bbeta <- (pvfm + 1) * exp(logfrailtypar)
    dist_id <- 2L
  }
  
  list(alpha = alpha, bbeta = bbeta, dist = dist_id)
}





















laplace_transform <- function(x, distribution) {
  # if(missing(.distribution) & missing())
  if(!inherits(distribution, "emfrail_dist"))
    stop("distribution argument misspecified; see ?emfrail_dist()")
  
  getpars <- dist_to_pars(distribution$dist, log(distribution$frailtypar), distribution$pvfm)
  
  if(getpars$dist == 0L) {
    L <- with(getpars, (bbeta / (bbeta + x))^alpha)
  }
  
  if(getpars$dist == 1L) {
    L <- with(getpars, exp(-1 * x^bbeta))
  }
  
  if(getpars$dist == 2L) {
    L <- with(getpars, exp(-alpha * sign(distribution$pvfm) * (1 - (bbeta / (bbeta + x))^distribution$pvfm )))
  }
  
  L
  
}


emfrail_pll <- function(formula, data,
                        distribution = emfrail_dist(),
                        values) {
  sapply(values, function(fp) {
    -emfrail(formula = formula,
             data = data,
             distribution =  emfrail_dist(dist = distribution$dist,
                                          theta = fp,
                                          pvfm = distribution$pvfm,
                                          left_truncation = distribution$left_truncation),
             control = emfrail_control(opt_fit = FALSE))
  })
  
  
  
}



