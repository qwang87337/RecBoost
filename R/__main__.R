library(dplyr)
library(survival)
library(Rcpp)
library(dplyr)
library(MASS)
library(mboost)

# import required functions
source("rCoxPH.R")
source("plloss.R")
source("em_arguments.R")
source("em_aux.R")
source("fast_Estep.R")
source("emboost_fit.R")

# read the sample data
df <- read.csv("sim_rec_dat.csv")

# in this sample, we have 100 non-informative covariates (V1,...V100) + 2 informative covariates (X1, X2)
num_covariates <- 100


start_time <- Sys.time()




# setups
unique_ids <- unique(df$id)
n.iteration <- 100
nu=0.1
distribution = emfrail_dist(theta = 2)
control = emfrail_control(se = F, em_control = list(eps = 0.001,
                                                    maxit = 20,
                                                    fast_fit = TRUE,
                                                    verbose = F,
                                                    upper_tol = exp(10),
                                                    lik_tol = 1),
                          nlm_control = list(stepmax = 1, gradtol = 1e-5, steptol = 1e-5,
                                             iterlim = 100
                          ))







#-----------------------------stability selection---------------------------#


all_covariates <- c("X1", "X2", paste0("V", 1:num_covariates))

# Initialize count vector
selection_count <- setNames(rep(0, length(all_covariates)), all_covariates)

# here we do resampling for stability selection 20 times for the purpose of demonstration. In 
# practice, it should be set on 100

n.resap <- 20
for (j in 1:n.resap) {
  
  sampled_ids <- sample(unique_ids, size = floor(0.75 * length(unique_ids)), replace = FALSE)
  
  # Subset the data to include only those subjects
  df_subsample <- df[df$id %in% sampled_ids, ]
  
  
  model = FALSE; model.matrix = FALSE
  sum_expression <- paste0("V", 1:num_covariates, collapse = " + ")
  fml <- as.formula(paste0("Surv(start, stop, status) ~ X1+X2+cluster(id)+", sum_expression))
  
  
  
  
  formula <- fml
  data=df_subsample
  nfrailty <- nrow(data)
  
  
  
  
  
  frail <- data[,c("id", "cen")]
  frail$ni <- data$num_rec
  frail$y <- Surv(data$start, data$stop, data$status)
  
  # check some prior setup
  if(isTRUE(control$em_control$fast_fit)) {
    if(!(distribution$dist %in% c("gamma", "pvf"))) {
      #message("fast_fit option only available for gamma and pvf with m=-1/2 distributions")
      control$em_control$fast_fit <- FALSE
    }
    
    # version 0.5.6, the IG fast fit gets super sensitive at small frailty variance...
    if(distribution$dist == "pvf")
      control$em_control$fast_fit <- FALSE
    
  }
  
  Call <- match.call()
  
  if(missing(formula) | missing(data)) stop("Missing arguments")
  
  cluster <- function(x) x
  terminal <- function(x) x
  strata <- function(x) x
  
  mf <- model.frame(formula, data)
  
  # Identify the positions of cluster and the ID column
  pos_cluster <- grep("cluster", names(mf))
  if(length(pos_cluster) != 1) stop("misspecified or non-specified cluster")
  id <- mf[[pos_cluster]]
  
  # find terminals(we assume no terminal so it can be deleted if want)
  pos_terminal <- grep("terminal", names(mf))
  if(length(pos_terminal) > 1) stop("misspecified terminal()")
  
  # find strata(we assume no strata so it can be deleted if want)
  pos_strata <- grep("strata", names(mf))
  if(length(pos_strata) > 0) {
    if(length(pos_strata) > 1) stop("only one strata() variable allowed")
    strats <- as.numeric(mf[[pos_strata]])
    label_strats <- levels(mf[[pos_strata]])
  } else {
    # else, everyone is in the same strata
    strats <- NULL
    label_strats <- "1"
  }
  
  
  
  
  
  Y <- mf[[1]]
  if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y) != 3) {              # col 1 the start time, col 2 the end time, col3 the outcome
    # making it all in (tstart, tstop) format
    Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
  }
  
  # create a matrix for explanatory variables
  X1 <- model.matrix(formula, data)
  
  pos_cluster_X1 <- grep("cluster", colnames(X1))
  pos_terminal_X1 <- grep("terminal", colnames(X1))
  pos_strata_X1 <- grep("strata", colnames(X1))
  
  
  
  # X is explanatory variable matrix (X1) without variables of cluster, terminal, strata and intercept
  X <- X1[,-c(1, pos_cluster_X1, pos_terminal_X1, pos_strata_X1), drop=FALSE]
  
  
  
  
  data2 <- cbind(Y,X)
  data2 <- as.data.frame(data2)
  
  colnames(data2)[1] <- "time1"
  colnames(data2)[2] <- "time2"
  colnames(data2)[3] <- "events"
  
  formula_str <- as.character(formula)
  formula_str <- gsub("\\+ cluster\\([^\\)]+\\)", "", formula_str)
  formula2 <- paste("Surv(time1, time2, events)",
                    formula_str[1],
                    formula_str[3])
  
  formula2 <- as.formula(formula2)
  
  
  
  # calculate weight for each observation
  freq <- table(id)
  weights <- 1 / freq[as.character(id)]
  weights <- as.numeric(weights)
  
  
  
  
  
  
  mb <-  glmboost(formula2, 
                  data = data2, family = rCoxPH(logz=rep(0, nfrailty), frail),
                  weights = weights,
                  control = boost_control(mstop = n.iteration))

  
  bnames <- names(unlist(mb$coef()))
  coef <- unlist(mb$coef())
  coef <- coef[names(coef) != "(Intercept)"]
  if ("(Intercept)" %in% bnames) {
    bnames <- bnames[bnames != "(Intercept)"]
  }
  X_sub <- X[,bnames]
  
  x2 <- matrix(rep(0, ncol(X)), nrow = 1, dimnames = list(123, dimnames(X)[[2]]))
  x2 <- scale(x2, center = apply(X, 2, mean), scale = FALSE)

  newrisk <- exp(c(x2[,bnames] %*% coef))
  exp_g_x <- exp(coef) %*% t(X_sub)                                   
  g <- coef                                                      
  g_x <- t(coef %*% t(X_sub))                                         
  explp <- exp(mb$predict()) # these are with centered covariates
  order_id <- match(id, unique(id))
  
  # number of event per cluster
  nev_id <- as.numeric(rowsum(Y[,3], order_id, reorder = FALSE))
  names(nev_id) <- unique(id)

  ord_tstop <- match(Y[,2], sort(unique(Y[,2])))     # ordered start time                                    
  ord_tstart <- match(Y[,1], sort(unique(Y[,1])))    # ordered stop time                                    
  
  nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
  esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))
  
  
  death <- (Y[, 3] == 1)
  nevent <- as.vector(rowsum(1 * death, Y[, ncol(Y) - 1]))
  time <- sort(unique(Y[,2])) 

  etime <- c(0, sort(unique(Y[, 1])),  max(Y[, 1]) + min(diff(time))/2)
  indx <- findInterval(time, etime, left.open = TRUE) # left.open  = TRUE is very important
  
  indx2 <- findInterval(Y[,1], time)
  
  time_to_stop <- match(Y[,2], time)
  
  atrisk <- list(death = death, nevent = nevent, nev_id = nev_id,
                 order_id = order_id,
                 time = time, indx = indx, indx2 = indx2,
                 time_to_stop = time_to_stop,
                 ord_tstart = ord_tstart, ord_tstop = ord_tstop,
                 strats = NULL)
  
  nrisk <- nrisk - c(esum, 0,0)[indx]
  if(newrisk == 0) warning("Hazard ratio very extreme; please check (and/or rescale) your data")
  haz <- nevent/nrisk * newrisk
  basehaz_line <- haz[atrisk$time_to_stop]
  cumhaz <- cumsum(haz)
  cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
  cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
  cumhaz_line <- (cumhaz[atrisk$time_to_stop] - c(0, cumhaz)[atrisk$indx2 + 1]) * explp / newrisk
  Cvec <- rowsum(cumhaz_line, order_id, reorder = FALSE)
  Cvec_lt <- 0 * Cvec           # Assuming no left truncation
  
  logz=rep(0, nfrailty)
  loglik <-  sum(plloss(Y, mb$predict(), frail, weights))

  
  outer_m <- do.call(nlm, args = c(list(f = emboost_fit,
                                        p = log(distribution$theta),
                                        hessian = TRUE,
                                        dist = distribution$dist,
                                        pvfm = distribution$pvfm,
                                        Y = Y, Xmat = X,
                                        atrisk = atrisk,
                                        basehaz_line = basehaz_line,
                                        mcox = list(coefficients = g, loglik = loglik),  # a "fake" cox model
                                        bnames = bnames,
                                        formulas = formula2,
                                        frail = frail, weights = weights,
                                        n.iteration = n.iteration,
                                        nu=nu,
                                        Cvec = Cvec,
                                        lt = distribution$left_truncation,
                                        Cvec_lt = Cvec_lt, se = FALSE,
                                        em_control = control$em_control),
                                   control$nlm_control))                                

  
  inner_m <- emboost_fit(logfrailtypar = outer_m$estimate, bnames = bnames,
                         dist = distribution$dist, pvfm = distribution$pvfm,
                         Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                         mcox = list(coefficients = g, loglik = mb$logLik()),  # a "fake" cox model
                         frail = frail, weights = weights,
                         Cvec = Cvec, lt = distribution$left_truncation,
                         Cvec_lt = Cvec_lt, se = FALSE,
                         formulas = formula2,
                         em_control = control$em_control,
                         n.iteration=n.iteration,
                         nu=nu,
                         return_loglik = FALSE)   
  

  selected_vars <- names(inner_m$coef)
  selection_count[selected_vars] <- selection_count[selected_vars] + 1

  
}


# get the variables that are selected 90% of the resampling.
selected_90plus <- names(selection_count[selection_count >= n.resap*0.9])

# Collapse into formula-like string
formula_sel <- paste(selected_90plus, collapse = " + ")


#--------------------------stability selection done--------------------------------#





model = FALSE; model.matrix = FALSE
# sum_expression <- paste0("V", 1:num_covariates, collapse = " + ")
fml <- as.formula(paste0("Surv(start, stop, status) ~ ", formula_sel, "+cluster(id)"))


formula <- fml
data=df
nfrailty <- nrow(data)


frail <- data[,c("id", "cen")]
frail$ni <- data$num_rec
frail$y <- Surv(data$start, data$stop, data$status)

# check some prior setup
if(isTRUE(control$em_control$fast_fit)) {
  if(!(distribution$dist %in% c("gamma", "pvf"))) {
    #message("fast_fit option only available for gamma and pvf with m=-1/2 distributions")
    control$em_control$fast_fit <- FALSE
  }
  
  # version 0.5.6, the IG fast fit gets super sensitive at small frailty variance...
  if(distribution$dist == "pvf")
    control$em_control$fast_fit <- FALSE
  
}

Call <- match.call()

if(missing(formula) | missing(data)) stop("Missing arguments")

cluster <- function(x) x
terminal <- function(x) x
strata <- function(x) x

mf <- model.frame(formula, data)

# Identify the positions of cluster and the ID column
pos_cluster <- grep("cluster", names(mf))
if(length(pos_cluster) != 1) stop("misspecified or non-specified cluster")
id <- mf[[pos_cluster]]

# find terminals(we assume no terminal so it can be deleted if want)
pos_terminal <- grep("terminal", names(mf))
if(length(pos_terminal) > 1) stop("misspecified terminal()")

# find strata(we assume no strata so it can be deleted if want)
pos_strata <- grep("strata", names(mf))
if(length(pos_strata) > 0) {
  if(length(pos_strata) > 1) stop("only one strata() variable allowed")
  strats <- as.numeric(mf[[pos_strata]])
  label_strats <- levels(mf[[pos_strata]])
} else {
  # else, everyone is in the same strata
  strats <- NULL
  label_strats <- "1"
}





# Y is basically the response variable in survival analysis, the output from Surv() function
Y <- mf[[1]]
if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
if(ncol(Y) != 3) {              # col 1 the start time, col 2 the end time, col3 the outcome
  # making it all in (tstart, tstop) format
  Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
}

# create a matrix for explanatory variables
X1 <- model.matrix(formula, data)

# no idea why we have this. Maybe cab delete it
pos_cluster_X1 <- grep("cluster", colnames(X1))
pos_terminal_X1 <- grep("terminal", colnames(X1))
pos_strata_X1 <- grep("strata", colnames(X1))



# X is explanatory variable matrix (X1) without variables of cluster, terminal, strata and intercept
X <- X1[,-c(1, pos_cluster_X1, pos_terminal_X1, pos_strata_X1), drop=FALSE]


# mcox the usual normal COX model using variables in X only, without cluster
# mcox <- survival::agreg.fit(x = X, y = Y, strata = strats, offset = NULL, init = NULL,
#                             control = survival::coxph.control(),
#                             weights = NULL, method = "breslow", rownames = NULL)
# mcox$coefficients

data2 <- cbind(Y,X)
data2 <- as.data.frame(data2)

colnames(data2)[1] <- "time1"
colnames(data2)[2] <- "time2"
colnames(data2)[3] <- "events"

formula_str <- as.character(formula)
formula_str <- gsub("\\+ cluster\\([^\\)]+\\)", "", formula_str)
formula2 <- paste("Surv(time1, time2, events)",
                  formula_str[1],
                  formula_str[3])

formula2 <- as.formula(formula2)



# calculate weight for each observation
freq <- table(id)
weights <- 1 / freq[as.character(id)]
weights <- as.numeric(weights)






mb <-  glmboost(formula2, 
                data = data2, family = rCoxPH(logz=rep(0, nfrailty), frail),
                weights = weights,
                control = boost_control(mstop = n.iteration))
summary(mb)

bnames <- names(unlist(mb$coef()))
coef <- unlist(mb$coef())
coef <- coef[names(coef) != "(Intercept)"]
if ("(Intercept)" %in% bnames) {
  bnames <- bnames[bnames != "(Intercept)"]
}
X_sub <- X[,bnames]



# x2 gives the mean from mcox
x2 <- matrix(rep(0, ncol(X)), nrow = 1, dimnames = list(123, dimnames(X)[[2]]))
x2 <- scale(x2, center = apply(X, 2, mean), scale = FALSE)
# x2 gives the mean from mcox
newrisk <- exp(c(x2[,bnames] %*% coef))
exp_g_x <- exp(coef) %*% t(X_sub)                                   
# the coefficients for variables from mcox without frailty                   
g <- coef                                                      
# g_x is a vector for the predicted hazard for each patient                  
g_x <- t(coef %*% t(X_sub))                                         
# explp the exponential of the predicted hazard for each obs
explp <- exp(mb$predict()) # these are with centered covariates
# match() returns the positions of the first match of id in unique(id)
order_id <- match(id, unique(id))
# number of event per cluster
nev_id <- as.numeric(rowsum(Y[,3], order_id, reorder = FALSE))
names(nev_id) <- unique(id)
# nrisk has the sum with every tstop and the sum of elp at risk at that tstop
# esum has the sum of elp who enter at every tstart
# indx groups which esum is right after each nrisk;
# the difference between the two is the sum of elp really at risk at that time point. 


# Assuming we have no strata

# we are using the following part to calculate cumulative hazard                      
ord_tstop <- match(Y[,2], sort(unique(Y[,2])))     # ordered start time                                    
ord_tstart <- match(Y[,1], sort(unique(Y[,1])))    # ordered stop time                                    


# the cumulative sum hazard
# Y[, ncol(Y) - 1] is the end time of obs.
# rowsum(explp, Y[, ncol(Y) - 1]): calculate the sum of explp, grouped by end time
# rev(...) flip what we got (1st to last, and so on...)
# cumsum() sums the things in () 
# so nrisk is a vector, with the sum of (predicted) hazard at each end time.
nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))


# Y[, 1] is the start time of each obs
# similar to nrisk, but it's the sum of predicted hazard at each start time
esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))

# logical indicator for death of each obs
death <- (Y[, 3] == 1)

# number of death at each time point
nevent <- as.vector(rowsum(1 * death, Y[, ncol(Y) - 1]))

#unique stop time
time <- sort(unique(Y[,2])) # unique tstops

# Y[,1] is the start time of each obs
# diff(x) gives a vector, resulting from x[2] - x[1], x[3] - x[2], so on...
# etime, a vector starting from 0, to each ordered start time, finally end with 
# a value = max(Y[, 1]) + min(diff(time))/2 (no idea why we need this)

# etime is a vector of unique "START" time with some small adjustment?
etime <- c(0, sort(unique(Y[, 1])),  max(Y[, 1]) + min(diff(time))/2)

# NOTE: length of indx is the length of unique end time from Y
# findInterval(x,y) categorize x into groups based on the interval you defined in y  
# try cbind(time, indx), which gives more intuition
indx <- findInterval(time, etime, left.open = TRUE) # left.open  = TRUE is very important


# nrisk has the sum with every tstop and the sum of elp at risk at that tstop
# esum has the sum of elp who enter at every tstart
# indx groups which esum is right after each nrisk;
# the difference between the two is the sum of elp really at risk at that time point.


# this gives for every tstart (line variable), after which event time did it come
indx2 <- findInterval(Y[,1], time)

# looks the same thing to the stop time of each obs. 
time_to_stop <- match(Y[,2], time)

atrisk <- list(death = death, nevent = nevent, nev_id = nev_id,
               order_id = order_id,
               time = time, indx = indx, indx2 = indx2,
               time_to_stop = time_to_stop,
               ord_tstart = ord_tstart, ord_tstop = ord_tstop,
               strats = NULL)

# nrisk is the sum of predicted hazard at each unique end time, 
# esum is the sum pf predicted hazard at each start time
# so this new nrisk is the difference of sum hazard between end time and start time? 
nrisk <- nrisk - c(esum, 0,0)[indx]
if(newrisk == 0) warning("Hazard ratio very extreme; please check (and/or rescale) your data")
haz <- nevent/nrisk * newrisk
basehaz_line <- haz[atrisk$time_to_stop]
cumhaz <- cumsum(haz)
cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
cumhaz_line <- (cumhaz[atrisk$time_to_stop] - c(0, cumhaz)[atrisk$indx2 + 1]) * explp / newrisk
Cvec <- rowsum(cumhaz_line, order_id, reorder = FALSE)
Cvec_lt <- 0 * Cvec           # Assuming no left truncation

logz=rep(0, nfrailty)
loglik <-  sum(plloss(Y, mb$predict(), frail, weights))





# loglik <-  sum((log(basehaz_line) + g_x)[Y[,3] == 1]) +
#   sum(Y[,3]) - sum((atrisk$nevent * log(atrisk$nevent))[atrisk$nevent > 0])



outer_m <- do.call(nlm, args = c(list(f = emboost_fit,
                                      p = log(distribution$theta),
                                      hessian = TRUE,
                                      dist = distribution$dist,
                                      pvfm = distribution$pvfm,
                                      Y = Y, Xmat = X,
                                      atrisk = atrisk,
                                      basehaz_line = basehaz_line,
                                      mcox = list(coefficients = g, loglik = loglik),  # a "fake" cox model
                                      bnames = bnames,
                                      formulas = formula2,
                                      frail = frail, weights = weights,
                                      n.iteration = n.iteration,
                                      nu=nu,
                                      Cvec = Cvec,
                                      lt = distribution$left_truncation,
                                      Cvec_lt = Cvec_lt, se = FALSE,
                                      em_control = control$em_control),
                                 control$nlm_control))                                


inner_m <- emboost_fit(logfrailtypar = outer_m$estimate, bnames = bnames,
                       dist = distribution$dist, pvfm = distribution$pvfm,
                       Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                       mcox = list(coefficients = g, loglik = mb$logLik()),  # a "fake" cox model
                       frail = frail, weights = weights,
                       Cvec = Cvec, lt = distribution$left_truncation,
                       Cvec_lt = Cvec_lt, se = FALSE,
                       formulas = formula2,
                       em_control = control$em_control,
                       n.iteration=n.iteration,
                       nu=nu,
                       return_loglik = FALSE)

res <- inner_m

# final results
print(res)
res$coef