plloss <- function(y, f, frail, w) {
  
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
    temp <- frail[frail$id==IDs[i],]                           # `temp` is the subset dataset of observation of subject k 
    
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