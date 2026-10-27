(c) Cambridge University Press 2007

-----------------------------------------------------
Problems and further results from Applied Asymptotics
-----------------------------------------------------

This file is AAproblems.R
+++++++++++++++++++++++++


## ----------------------------------------------------------------------
## Problem 40  
## ----------------------------------------------------------------------
library(nlreg)
#
log.time <- log(daphnia$time)
x0.seq <- seq(0, 0.2, 0.01)
logLik.x0 <- matrix(0, nrow = length(x0.seq), ncol = 2)
for( idx in seq(along = x0.seq) )
{
  x0 <- x0.seq[idx]
  conc.x0 <- ifelse(daphnia$conc - x0 > 0, daphnia$conc - x0, 0)
  fitted.mod <- 
    nlreg( log.time ~ log(b1+(b2-b1)/(1+(conc.x0/b4)^(2*b3))), 
           start = list(b1 = 0.05, b2 = 50, b3 = 3, b4 = 0.1), 
           hoa = TRUE )
  logLik.x0[idx,] <- c( x0, logLik(fitted.mod) )
}
plot( logLik.x0, type = 'l', ylab = 'profile log likelihood',
      xlab = expression(paste('threshold ', x[0])) )
## ======================================================================
