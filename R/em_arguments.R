emfrail_control <- function(opt_fit = TRUE,
                            se = TRUE,
                            se_adj = TRUE,
                            ca_test = TRUE,
                            lik_ci = TRUE,
                            lik_interval = exp(c(-3, 20)),
                            lik_interval_stable = exp(c(0, 20)),
                            nlm_control = list(stepmax = 1),
                            zph = FALSE,
                            zph_transform = "km",
                            em_control = list(eps = 0.0001,
                                              maxit = Inf,
                                              fast_fit = TRUE,
                                              verbose = FALSE,
                                              upper_tol = exp(10),
                                              lik_tol = 1)
) {
  # calculate SE as well
  
  # Here some checks
  
  if(isTRUE(lik_ci)) {
    if(length(lik_interval) != 2 | length(lik_interval_stable) != 2)
      stop("lik_interval must be of length 2")
    if(lik_interval[1] < exp(-7) | lik_interval[2] > exp(20))
      warning("extreme values for lik_interval, there might be some numerical trouble")
    # if(lik_ci_intervals$interval[2] != em_control$upper_tol)
    # message("it is good practice right hand side of the interval should be equal to upper_tol")
  }
  
  # make sure the defaults of these function are the same as those from the input!
  inner_c <- function(eps = 0.0001,
                      maxit = Inf,
                      fast_fit = TRUE,
                      verbose = FALSE,
                      upper_tol = exp(20),
                      lik_tol = 1) {
    list(eps = eps,
         maxit = maxit,
         fast_fit = fast_fit,
         verbose = verbose,
         upper_tol = upper_tol,
         lik_tol = lik_tol)
  }
  
  em_control <- do.call(inner_c, em_control)
  
  res <- list(opt_fit = opt_fit,
              se = se,
              se_adj = se_adj,
              ca_test = ca_test,
              lik_ci = lik_ci,
              lik_interval = lik_interval,
              lik_interval_stable = lik_interval_stable,
              zph = zph,
              zph_transform = zph_transform,
              nlm_control = nlm_control,
              em_control = em_control)
  attr(res, "class") <- c("emfrail_control")
  res
}




#' Distribution parameters for emfrail
#'
#' @param dist One of 'gamma', 'stable' or 'pvf'.
#' @param theta A starting value for the 'outer' maximization with respect to the frailty parameter \eqn{\theta}. Must be >0.
#' @param pvfm Only relevant if \code{dist = 'pvf'} is used. It determines which PVF distribution should be used. Must be  larger than -1 and not equal to 0.
#' @param left_truncation Logical. Whether the data set represents left truncated survival times.
#' @param basehaz A character string which determines how the baseline hazard is calculated. The default is "breslow", but other possible options are "weibull", "exponential" "gaussian", "logistic", "lognormal" or "loglogistic".
#' @return An object of the type \code{emfrail_dist}, which is mostly used to denote the
#' supported frailty distributions in a consistent way.
#' @export
#'
#' @details The \code{theta} argument must be positive. In the case of gamma or PVF, this is the inverse of
#'  the frailty variance, i.e. the larger the \code{theta} is,
#'  the closer the model is to a Cox model. When \code{dist = "pvf"} and \code{pvfm = -0.5}, the inverse Gaussian
#'  distribution is obtained. For the positive stable distribution, the \eqn{\gamma} parameter of the Laplace transform is
#'  \eqn{\theta / (1 + \theta)}, with the \eqn{alpha} parameter fixed to 1.
#'
#' @seealso \code{\link{emfrail}, \link{emfrail_control}}
#' @examples
#' emfrail_dist()
#' # Compound Poisson distribution:
#' emfrail_dist(dist = 'pvf', theta = 1.5, pvfm = 0.5)
#' # Inverse Gaussian distribution:
#' emfrail_dist(dist = 'pvf')

emfrail_dist <- function(dist = "gamma", theta = 2, pvfm = -1/2, left_truncation = FALSE, basehaz = "breslow") {
  
  if (!(dist %in% c("gamma", "stable", "pvf")))
    stop("frailty distribution must be one of gamma, stable, pvf")
  if (length(theta) != 1)
    stop("specify exactly 1 parameter (theta>0) for the frailty")
  if (theta <= 0)
    stop("frailty parameter (theta) must be positive")
  if (dist == "pvf" & (pvfm < -1 | pvfm == 0))
    stop("pvfm must be >-1 and not equal to 0")
  
  if(!is.logical(left_truncation)) stop("left_truncation must be TRUE or FALSE")
  
  if(basehaz != "breslow") {
    # here some check for this stuff so that it matches something
  }
  
  res <- list(dist = dist, theta = theta, pvfm = pvfm, left_truncation = left_truncation, basehaz = basehaz)
  attr(res, "class") <- c("emfrail_dist")
  return(res)
}