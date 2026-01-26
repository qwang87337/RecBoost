rCoxPH <- function(logz, frail) {
  
  plloss <- function(y, f, w) {
    
    w <- NULL # no need for weights in this
    
    frail$logz <- logz
    frail2 <- frail[!duplicated(frail$id), ]
    
    n <- length(unique(frail$id)) # number of subject
    IDs <- unique(frail$id)
    expf <- exp(f + frail$logz)
    expf <- aggregate(expf, by = list(frail$id), FUN = function(x) x[1])
    expf <- as.numeric(expf[,2])
    
    f2 <- aggregate(f, by = list(frail$id), FUN = function(x) x[1])
    f2 <- as.numeric(f2[,2])
    
    res <- numeric(n)
    
    t2 <- 0
    for (i in 1:n) {
      
      t2 <- 0
      ni <- as.numeric(frail2[frail2$id == IDs[i], ]$ni)         # find the number of event for subject i
      temp <- frail[frail$id==i,]                           # `temp` is the subset dataset of observation of subject k 
      
      if (ni == 0) {
        
        t2 <- t2 + 0
        
      } else {
        
        for (j in 1:ni) {
          
          indicator <- as.numeric(frail2$cen >= temp$y[,2][j]) # a logical vector telling you that out of 100 subjects, who are in the risk set at time t_kj
          term2 <- sum(indicator*expf)          # the denominator: summation of all predictions in the risk set 
          
          t2 <- t2 + f2[i] - log(term2)
        } 
      }
      res[i] <- t2
    }
    return(res)
  }
  
  Family(
    
    ngradient <- function(y, f, w = NULL) {
      w <- NULL # no need for weights in this
      
      frail$logz <- logz
      frail2 <- frail[!duplicated(frail$id), ]
      
      
      n <- length(unique(frail$id)) # number of subject
      
      f <- as.vector(f)
      f <- f + logz
      expf <- exp(f)
      
      
      
      frail$logz <- logz
      frail$expf <- expf
      frail2 <- frail[!duplicated(frail$id), ]
      frail3 <- frail[frail$y[,3]==1,]
      
      
      
      risk_matrix <- as.matrix(outer(frail2$cen, frail3$y[, 2], ">="))
      trisk <- t(risk_matrix)
      
      # this is a vector of denominator (after summation) at each event time (so length equals number of events)
      # this way of calculating denominator is 100% correct! trust yourself!!!!!
      denominator <- as.vector(trisk %*% frail2$expf)
      
      u = rep(0,n)
      n_events <- as.vector(frail2$ni)
      for (i in 1:n) {
        
        # a vector of length of number of events, indicating if the i-th subject's censoring time >= event time t_kj 
        numinator <- as.vector(risk_matrix[i,]*frail2$expf[i])
        
        u[i] <- n_events[i] - sum(numinator/denominator)
        
      }
      row_counts <- frail %>%
        group_by(id) %>%
        summarise(row_count = n())
      
      u <- rep(u, as.numeric(row_counts$row_count))
      return(u)
    },
    risk = risk <- function(y, f, w = 1) -sum(plloss(y, f, w), na.rm = TRUE),
    offset = function(y, w = 1) 0, ## perhaps use something different
    
    ## Note: offset cannot be computed from Cox Partial LH as
    ## PLH doesn't depend on constant
    name = "frailty Cox frailty Partial Likelihood for calendar time"
  )
}
