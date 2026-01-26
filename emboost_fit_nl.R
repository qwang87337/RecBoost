emboost_fit_nl <- function(logfrailtypar, dist, pvfm,
                        Y, Xmat, 
                        atrisk, 
                        basehaz_line, 
                        mcox = list(),
                        lp, frail, weights,
                        Cvec, lt = FALSE, Cvec_lt, # we need to start somewhere with the Cvec (E step comes first)
                        em_control, se,
                        return_loglik = TRUE,
                        formulas,
                        n.iteration) {
  

  pars <- dist_to_pars(dist, logfrailtypar, pvfm)
  if (isTRUE(em_control$verbose)) {
    print(paste0(#"dist=", pars$dist,
      "logfrailtypar= ", logfrailtypar,
      " / alpha=", pars$alpha,
      " / bbeta=", pars$bbeta))
  }
  
  if(logfrailtypar < -100) warning("theta virtually 0; try another starting value")
  
  g_x <- lp
  
  
  # if the logfrailtypar is large (i.e. frailty variance 0) then just return the Cox likelihood
  if(logfrailtypar > log(em_control$upper_tol)) {
    
    #message("Frailty parameter very large, frailty variance close to 0")
    loglik <- mcox$loglik[length(mcox$loglik)]
    
    
    if(isTRUE(return_loglik)) {
      if(isTRUE(em_control$verbose)) print(paste("loglik = ",loglik))
      return(-loglik)
    }
    
  }
  
  
  loglik_old = -Inf
  ncycles <- 0
  convergence <- FALSE
  
  # EM algorithm iteration part
  while(!isTRUE(convergence)) {
    
    if(isTRUE(em_control$fast_fit)) {
      e_step_val <- fast_Estep(Cvec,
                               Cvec_lt,
                               atrisk$nev_id,
                               alpha = pars$alpha,
                               bbeta = pars$bbeta,
                               pvfm = pvfm,
                               dist = pars$dist)
    } else {
      e_step_val <- Estep(Cvec,
                          Cvec_lt,
                          atrisk$nev_id,
                          alpha = pars$alpha,
                          bbeta = pars$bbeta,
                          pvfm = pvfm,
                          dist = pars$dist)
    }
    
    
    
    # logz is the result from log((e_step_val[,1] / e_step_val[,2]) given to each obs.
    # the obs from the same ID (cluster) has the same logz
    # based on the equation, I believe logz is the log of z, the random effect
    logz <- log((e_step_val[,1] / e_step_val[,2])[atrisk$order_id])
    
    frail$logz <- logz
    
    # assume no strata
    # this loglik is the log of likelihood given by equation (3) in
    # page 5 of Balan et al., 2019
    loglik <-  sum((log(basehaz_line) + g_x)[Y[,3] == 1]) + sum(e_step_val[,3]) +
      sum(Y[,3]) - sum((atrisk$nevent * log(atrisk$nevent))[atrisk$nevent > 0])
    
    mcox$loglik <- loglik
    
    # if this happens, then something is going very wrong
    if(loglik < loglik_old - em_control$lik_tol)
      warning(paste0("likelihood decrease of ", loglik - loglik_old ))
    
    
    print(loglik)
    # here is the break part of iteration
    if(abs(loglik - loglik_old) < em_control$eps) break
    
    loglik_old <- loglik
    
    
    # replace agreg.fit() with glmboost()
    # mcox <- survival::agreg.fit(x = Xmat, y = Y, strata = atrisk$strats, offset = logz, init = NULL,
    #                             control = survival::coxph.control(), weights = NULL,
    #                             method = "breslow", rownames = NULL)
    # 
    
    # mcox = list(coefficients = g, loglik = mcox$loglik)
    # seems I have to put Xmat and Y together to form a survival data...
    
    newdata <- cbind(Y[,1],Y[,2],Y[,3], Xmat)
    newdata <- as.data.frame(newdata)
    
    colnames(newdata)[1] <- "time1"
    colnames(newdata)[2] <- "time2"
    colnames(newdata)[3] <- "events"
    
    
    # formula_str1 <- "Surv(time1, time2, events)~"
    # formula_str2 <- paste(colnames(newdata)[-c(1,2,3)], collapse = "+")
    # new.form <- paste(formula_str1, formula_str2)
    # new.form <- as.formula(new.form)
    
    
    # need to add some code here that identify the type of 
    mboost_res <- gamboost(Surv(time1, time2, events) ~ bols(X1) + bbs(X2), 
                           data = newdata, 
                           family = rCoxPH(logz <- logz, frail = frail), 
                           weights = weights,
                           control = boost_control(mstop = n.iteration))
    
    
    
    ncycles <- ncycles + 1
    
    
    
    bnames <- names(unlist(mb$coef()))
    coef <- unlist(mb$coef())
    coef <- coef[names(coef) != "(Intercept)"]
    
    if ("(Intercept)" %in% bnames) {
      bnames <- bnames[bnames != "(Intercept)"]
    }
    
    
    
    
    
    pred <- mboost_res$predict()
    
    mcox = list(coefficients = coef, loglik = mboost_res$logLik())
    
    # g_x <- t(coef %*% t(Xmat[,bnames]))
    
    lp <- pred + as.numeric(mean(mb$predict()))
    g_x <- lp
    
    explp <- exp(lp)
    
    nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))
    nrisk <- nrisk - c(esum, 0,0)[atrisk$indx]
    
    haz <- atrisk$nevent/nrisk
    basehaz_line <- haz[atrisk$time_to_stop]
    cumhaz <- cumsum(haz)
    
    cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
    
    cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
    
    cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)
    
    
    # assuming not left truncation
    Cvec_lt <- 0 * Cvec
    
    
    Cvec <- rowsum(cumhaz_line * exp(lp), atrisk$order_id, reorder = FALSE)
    
    
    if(ncycles > em_control$maxit) {
      warning(paste("did not converge in ", em_control$maxit," iterations." ))
      break
    }
  }
  
  
  if(isTRUE(return_loglik)) {
    if(isTRUE(em_control$verbose)) print(paste("loglik = ",loglik))
    return(-loglik)
  }
  
  
  
  if(!isTRUE(return_loglik)) {
    
    res = list(
      model = mboost_res,
      coef = unlist(mboost_res$coef()),
      loglik = mboost_res$logLik(),
      logfrailtypar = logfrailtypar,
      sel.freq = summary(mboost_res)$selprob,
      hazard = haz,
      cumhaz_line = cumhaz_line,
      frail = e_step_val[,1] / e_step_val[,2],
      logz = logz
    )
    
    return(res)
    
  }
  
  
  
}
