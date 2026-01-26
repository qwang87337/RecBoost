emboost_fit <- function(logfrailtypar, dist, pvfm,
                            Y, Xmat, # id,  # this is some data stuff
                            atrisk, # a list with stuff that won't change in the EM
                            basehaz_line ,  # need for log-likelihood
                            mcox = list(),
                            bnames, frail, weights,
                            Cvec, lt = FALSE, Cvec_lt, # we need to start somewhere with the Cvec (E step comes first)
                            em_control, se,
                            return_loglik = TRUE,
                            formulas,
                            n.iteration, 
                            nu) {
  
  Xmat_sub <- Xmat[,bnames]
  pars <- dist_to_pars(dist, logfrailtypar, pvfm)
  if (isTRUE(em_control$verbose)) {
    print(paste0(#"dist=", pars$dist,
      "logfrailtypar= ", logfrailtypar,
      " / alpha=", pars$alpha,
      " / bbeta=", pars$bbeta))
  }
  
  if(logfrailtypar < -100) warning("theta virtually 0; try another starting value")
  
  g_x <- t(mcox$coefficients %*% t(Xmat_sub))
  lp <- g_x
  
  # if the logfrailtypar is large (i.e. frailty variance 0) then just return the Cox likelihood
  if(logfrailtypar > log(em_control$upper_tol)) {
    
    #message("Frailty parameter very large, frailty variance close to 0")
    loglik <- mcox$loglik[length(mcox$loglik)]
    
    
    if(isTRUE(return_loglik)) {
      if(isTRUE(em_control$verbose)) print(paste("loglik = ",loglik))
      return(loglik)
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
    
    
    
    newdata <- cbind(Y[,1],Y[,2],Y[,3], Xmat)
    newdata <- as.data.frame(newdata)
    
    colnames(newdata)[1] <- "time1"
    colnames(newdata)[2] <- "time2"
    colnames(newdata)[3] <- "events"
  

    # loglik <-  sum(plloss(y, lp, w))
    
    
    
    
    loglik <-  sum((log(basehaz_line) + g_x)[Y[,3] == 1]) + sum(e_step_val[,3]) +
      sum(Y[,3]) - sum((atrisk$nevent * log(atrisk$nevent))[atrisk$nevent > 0])
    
    

    
    mcox$loglik <- loglik
    
    # if this happens, then something is going very wrong
    if(loglik < loglik_old - em_control$lik_tol)
      warning(paste0("likelihood decrease of ", loglik - loglik_old ))
    
    
    print(paste0("-loglik=",-loglik))
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
    

    
    
    # formula_str1 <- "Surv(time1, time2, events)~"
    # formula_str2 <- paste(colnames(newdata)[-c(1,2,3)], collapse = "+")
    # new.form <- paste(formula_str1, formula_str2)
    # new.form <- as.formula(new.form)
    
    
   
    mboost_res <- glmboost(formulas, data = newdata, 
                           family = rCoxPH(logz <- logz, frail = frail), 
                           weights = weights,
                           control = boost_control(mstop = n.iteration, nu = nu))
    

    
    ncycles <- ncycles + 1
    
    
  
    bnames <- names(unlist(mboost_res$coef()))
    coef <- unlist(mboost_res$coef())
    coef <- coef[names(coef) != "(Intercept)"]
    
    if ("(Intercept)" %in% bnames) {
      bnames <- bnames[bnames != "(Intercept)"]
    }
    
    
    
    
    
    pred <- mboost_res$predict()
    
    mcox = list(coefficients = coef, loglik = mboost_res$logLik())
    
    g_x <- t(coef %*% t(Xmat[,bnames]))
    
    ncolumns_temp <- ncol(as.matrix(Xmat[,bnames]))
    
    if(ncolumns_temp == 1) {
      
      rescale0525 <- coef * mean(Xmat[,bnames])
      
    } else {rescale0525 <- as.numeric(coef %*% colMeans(Xmat[,bnames]))}
    
    
    # lp <- pred + as.numeric(coef %*% colMeans(Xmat[,bnames]))

    lp <- pred + rescale0525
    
    
    explp <- exp(lp)
    
    nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))
    nrisk <- nrisk - c(esum, 0,0)[atrisk$indx]
    
    haz <- atrisk$nevent/nrisk
    basehaz_line <- haz[atrisk$time_to_stop]
    cumhaz <- cumsum(haz)
    
    # cumhaz_0_line is cumhaz re-ordered by time...? 
    cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
    
    cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
    
    cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)
    
    
    # assuming not left truncation
    Cvec_lt <- 0 * Cvec
    
    
    Cvec <- rowsum(cumhaz_line * exp(g_x), atrisk$order_id, reorder = FALSE)
    
    
    if(ncycles > em_control$maxit) {
      warning(paste("did not converge in ", em_control$maxit," iterations." ))
      break
    }
  }
  
  
  if(isTRUE(return_loglik)) {
    if(isTRUE(em_control$verbose)) print(paste("-loglik = ",-loglik))
    return(-loglik)
  }
  

  
  if(!isTRUE(return_loglik)) {
    
    res = list(
      
      model = mboost_res,
      coef = unlist(mboost_res$coef()),
      loglik = loglik,
      logfrailtypar = logfrailtypar,
      sel.freq = summary(mboost_res)$selprob,
      cumhaz = cumhaz,
      cumhaz_0_line <- cumhaz_0_line,
      cumhaz_line = cumhaz_line,
      frail = e_step_val[,1] / e_step_val[,2],
      logz = logz
    )
    
    return(res)
  
  }


  
}




