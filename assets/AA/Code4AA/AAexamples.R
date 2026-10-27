(c) Cambridge University Press 2007

--------------------------------------------------
Examples and case studies from Applied Asymptotics
--------------------------------------------------

This file is AAexamples.R
+++++++++++++++++++++++++

 
===============================
CHAPTER 3: Simple Illustrations
=============================== 


## ----------------------------------------------------------------------
## Code 3.1: Code to compare first order and third order approximations 
##           for a Cauchy observation.
## ----------------------------------------------------------------------
score <- function(y) sqrt(2) * (y-theta)/(1+(y-theta)^2)
root <- function(y) sign(y-theta) * sqrt( 2*(log(1+(y-theta)^2)) )
wald <- function(y) sqrt(2) * (y-theta)
lugr <- function(y) pnorm( root(y) ) + 
                      dnorm( root(y) ) * ( 1/root(y) - 1/score(y) )
rstar <- function(y) pnorm( root(y) + log(score(y)/root(y))/root(y) )

y <- seq(from = 0, to = 15, by = 0.1)
theta <- 0

plot(y, pnorm(wald(y)), xlab = "y", ylab = "Distribution function", 
     type = "l", lty = 6)
lines(y, pnorm(score(y)), lty = 5)
lines(y, pnorm(root(y)), lty = 4)
lines(y, lugr(y), lty = 3)
lines(y, rstar(y), lty = 2)
lines(y, pcauchy(y), lty = 1)
legend(9, 0.9, c("Exact", "rstar", "Lugannani-Rice", "root", "score", 
                 "Wald"), lty = c(1,2,3,4,5,6), bty="n")
## ======================================================================


## ----------------------------------------------------------------------
## Code 3.2: R code to compute exact and third order approximate 
##           significance levels for a Poisson observation.
## ----------------------------------------------------------------------
root <- function(th, y) { sign(y-th) * sqrt( 2*(y*log(y/th)-y+th) ) }
wald <- function(th, y) { sqrt(y)*log(y/th) }
lugr <- function(th, y = 27) {
  pnorm( root(th, y) ) + 
    dnorm( root(th, y) ) * ( 1/root(th, y) - 1/wald(th, y) ) }
rstar <- function(th, y = 27) {
  pnorm( root(th, y) + log( wald(th, y)/root(th, y) )/root(th, y) ) }

b <- 6.7
mu <- seq(from = -2, to = 40, by = 0.01)
y0 <- 27

# Create Figure 3.3

plot(mu, ppois(y0, b+mu, lower.tail = F),
     type = "n", xlab = "mu", ylab = "Significance function")
lines(mu, ppois(y0-1, b+mu, lower.tail = F), lty = 2)
lines(mu, ppois(y0, b+mu, lower.tail = F), lty = 2)
lines(mu, ppois(y0, b+mu, lower.tail = F) + 0.5 * dpois(y0, b+mu), 
      lty = 3)
lines(mu, 1 - rstar(b+mu))

# Compute significance levels

ppois( y0-1, b, lower.tail = F )
1 - rstar(b, y0 - 0.5)
ppois( y0, b, lower.tail = F )
1 - rstar(b, y0 + 0.5)
ppois( y0, b, lower.tail = F ) + 0.5 * dpois(y0, b)
1 - rstar(b)
pnorm( (y0-1-b)/sqrt(b), lower.tail = F )
pnorm( (y0-b)/sqrt(b), lower.tail = F )
pnorm( (y0-1-b)/sqrt(y0), lower.tail = F ) 
pnorm( (y0-b)/sqrt(y0), lower.tail = F )

# Create Figure 3.2

y <- seq(from = -1, to = 30, by = 1)
plot(y, ppois(y, b), type = "n", xlim = c(0, 30), xlab = "y", 
     ylab = "Cumulative distribution function")
abline(v = y + 1.5, col = "gray")
lines(y, ppois(y, b), type = "s")
x <- seq(from = 0, to = 30, length = 1000)
lines(x, rstar(b, x), lty = 2)
lines(x, rstar(b, x+0.5), lty = 3)
## ======================================================================


## ----------------------------------------------------------------------
## Code 3.3: Approximate conditional inference for the Hubble data.
## ----------------------------------------------------------------------
hubble <- data.frame( r = c(1, 5), m = c(19, 7), x = c(0, 1) )
hubble

summary( hubble.fit <- glm( cbind(r, m-r) ~ x, binomial, data = hubble ) )

library( cond )
summary( cond.glm( hubble.fit, offset = x ), test = 0 )
## ======================================================================


## ----------------------------------------------------------------------
## Code 3.4: Higher order inference for the ratio for the Hubble data.
## ----------------------------------------------------------------------
nlogL <- function(psi, lam, d) 
{  
  m <- min(1, exp(-psi))
  p0 <- m / (1 + exp(-lam))
  p1 <- p0 * exp(psi)
  - sum( dbinom(d$r, size = d$m, prob = c(p0, p1), log = T) )
}

make.V <- function(th, d)  NULL   # V is not needed in this case

phi <- function(th, V, d)
{  
  psi <- th[1]
  lam <- th[-1]
  m <- min(1, exp(-psi))
  p0 <- m / (1 + exp(-lam))
  p1 <- p0 * exp(psi)
  phi1 <- log( p0/(1-p0) )
  phi2 <- log( p1/(1-p1) )
  c(phi1, phi2)
}

hubble <- data.frame( r = c(1, 5), m = c(19, 7), x = c(0, 1) )

library( fraser.reid )
( hubble.fr <- fraser.reid( 0, nlogL, phi, make.V, 
                            th.init = c(0, 0.1), hubble ) ) 
## ======================================================================


## ----------------------------------------------------------------------
## Code 3.5: Higher order inference for the log ratio of means under the
##           log-normal model for the cost data.
## ----------------------------------------------------------------------
# Make data frame and basic functions

cost <- data.frame(
          f = factor( c(rep(1, 13), rep(2, 18)) ),
          y = c( 30,172,210,212,335,489,651,1263,1294,1875,2213,
                 2998,4935,121,172,201,214,228,261,278,279,351,561,
                 622,694,848,853,1086,1110,1243,2543 ) )
  
nlogL <- function(psi, lam, d) {
#  d = data
  s1 <- exp(lam[2])
  m2 <- lam[1]
  s2 <- exp(lam[3])
  m1 <- psi + m2 + s2^2/2 - s1^2/2
  - sum( dnorm( log(d$y), mean = ifelse(d$f==1, m1, m2),
                sd = ifelse(d$f==1, s1, s2), log = T ) )
}

phi <- function(th, V, d) {
  psi <- th[1]
  lam <- th[-1]
  s1 <- exp(lam[2])
  m2 <- lam[1]
  s2 <- exp(lam[3])
  m1 <- psi + m2 + s2^2/2 - s1^2/2
  c( m1/s1^2, 1/s1^2, m2/s2^2, 1/s2^2 ) 
}

make.V <- function(th, d)  NULL

library( fraser.reid )
cost.lnorm.rat <- fraser.reid( psi = NULL, nlogL = nlogL, phi = phi, 
                               make.V = make.V, 
                               th.init = c(0, 5, 2, 5), d = cost )
  
plot.fr( cost.lnorm.rat )
lik.ci( cost.lnorm.rat )
## ======================================================================


## ----------------------------------------------------------------------
## Code 3.6: Bootstrap inference for the difference and ratio for the 
##           cost data.
## ----------------------------------------------------------------------
# Basic bootstrapping functions
  
stat <- function(data) {  
# computation of difference and ratio of means, and estimated variances
  m <- tapply( data$y, data$f, mean )
  v <- tapply( data$y, data$f, var )
  n <- tapply( data$y, data$f, length )
  c( m[1] - m[2], v[1]/n[1] + v[2]/n[2], log( m[1]/m[2] ), 
     v[1]/(m[1]^2 * n[1]) + v[2]/(m[2]^2 * n[2]) ) 
}

cost <- data.frame(
          f = factor( c(rep(1, 13), rep(2, 18)) ),
          y = c( 30,172,210,212,335,489,651,1263,1294,1875,2213,
                 2998,4935,121,172,201,214,228,261,278,279,351,561,
                 622,694,848,853,1086,1110,1243,2543 ) )

cost.f <- function(data, i = 1:nrow(data))  stat(cost[i,])  

library( boot )

# Nonparametric bootstrap

boot.ci( cost.boot <- boot( cost, cost.f, R = 10000, 
                            strata = cost$f ) )  # for difference
plot(cost.boot)  
boot.ci( cost.boot, index = c(3,4) )  # results for log ratio of means
plot(cost.boot, index = 3)  

# Parametric bootstrap using exponential model, giving results for 
# difference of means

cost.mle <- c( tapply(cost$y, cost$f, mean), 
               tapply(cost$y, cost$f, length) )
cost.sim <- function(data, mle) {   
  d <- data
  d$y <- c( rexp(mle[3], 1/mle[1]), rexp(mle[4], 1/mle[2]) )
  d 
}
boot.ci( 
  cost.para <- boot( cost, stat, R = 10000, ran.gen = cost.sim, 
                     mle = cost.mle, sim = "parametric" ), 
        type = c("norm", "stud") )

# Parametric bootstrap using log-normal model, giving results for 
# difference of means

v <- tapply( log(cost$y), cost$f, var )
n <- tapply( cost$y, cost$f, length )
sig <- (n-1) * v/n
cost.mle.lnorm <- c( tapply(log(cost$y), cost$f, mean), sqrt(sig), n )
cost.sim.lnorm <- function(data, mle) {
  d <- data
  d$y <- exp( c( rnorm(mle[5], mean = mle[1], sd = mle[3]),
                 rnorm(mle[6], mean = mle[2], sd = mle[4]) ) )
  d 
}
boot.ci( 
  cost.para.lnorm <- boot( cost, stat, R = 10000,
                           ran.gen = cost.sim.lnorm,
                           mle = cost.mle.lnorm, sim = "parametric" ), 
        type = c("norm", "stud") )
## ======================================================================


========================
CHAPTER 4: Discrete Data
========================


## ----------------------------------------------------------------------
## Code 4.1: Approximate conditional inference for the urine data.
## ----------------------------------------------------------------------
library( cond )

urine.glm <- 
  glm( formula = r ~ I(100*(gravity-1))+ph+osmo+conduct+urea+calc, 
       family = binomial, data = urine ) 
summary( urine.glm )

urine.cond.urea <- cond( urine.glm, offset = urea )
coef( urine.cond.urea )
summary( urine.cond.urea, coef = FALSE )

urine.cond.calc <- cond( urine.glm, offset = calc )
coef( urine.cond.calc )
summary( urine.cond.calc, coef = FALSE )
## ======================================================================


## ----------------------------------------------------------------------
## Code 4.2: R code for analysis of risk difference for the multiple 
##           myeloma data.
## ----------------------------------------------------------------------
n1 <- c(4, 4,2,2,2,3,2,5,2,2,3,2,4,3,4,12,2,3,4,3,4)
y1 <- c(3, 3,2,2,2,1,2,1,2,0,3,2,1,2,2, 4,1,3,1,0,2)
n2 <- c(3,11,3,2,3,3,3,4,3,3,3,2,5,4,6, 9,3,4,3,2,5)
y2 <- c(1, 8,2,2,0,2,2,4,2,2,3,0,1,2,4, 3,2,1,2,0,1)

risk <- data.frame( n1, y1, n2, y2 )

nlogL <- function(psi, lam, d) 
{
  xi0 <- max(0, -psi)
  xi1 <- min(1, 1 - psi)
  xi <- xi0 + (xi1-xi0) / (1 + exp(-lam))
  p1 <- psi + xi
  p2 <- xi
  - sum( d$y1*log(p1) + (d$n1-d$y1) * log(1-p1) + 
         d$y2*log(p2) + (d$n2-d$y2) * log(1-p2) )
}

make.V <- function(th, d)  NULL

phi <- function(th, V, d)
{
  n <- nrow(d)
  out <- rep(0, n + 1)
  psi <- th[1]
  lam <- th[-1]
  xi0 <- max(0, -psi)
  xi1 <- min(1, 1 - psi)
  xi <- xi0 + (xi1-xi0) / (1 + exp(-lam))
  p1 <- psi + xi
  p2 <- xi
  out1 <- d$n1 * log( p1/(1-p1) )
  out2 <- d$n2 * log( p2/(1-p2) )
  c( sum(out1), apply( cbind(out1, out2), 1, sum ) )
}

bad <- c(4, 11, 20)
psi.val <- seq(from = -0.2, to = 0.3, length = 50)
risk.init <- rep(0, 22) 

library( fraser.reid )
risk.fr <- 
  fraser.reid( psi = psi.val, nlogL = nlogL, phi = phi, make.V = make.V, 
               th.init = risk.init[-bad], d = risk[-bad,] )
## ======================================================================


===============================================
CHAPTER 5: Regression with Continuous Responses
===============================================


## ----------------------------------------------------------------------
## Code 5.1: R code for analysis of nuclear power station data.
## ----------------------------------------------------------------------
library( marg )

# Fit normal-theory linear model and examine its contents

nuc.norm <- lm( log(cost) ~ date + log(cap) + NE + CT + log(N) + PT, 
                data = nuclear )
summary( nuc.norm )

# Fit linear model with t errors and 4 df and examine its contents

nuc.t4 <- rsm( log(cost) ~ date + log(cap) + NE + CT + log(N) + PT, 
               data = nuclear, family = student(4) )
summary( nuc.t4 )
plot( nuc.t4 )

# Conditional analysis for partial turnkey guarantee

nuc.t4.pt <- cond( nuc.t4, offset = PT )
summary( nuc.t4.pt )
plot( nuc.t4.pt )

# For conditional analysis for other covariates, replace PT by 
# log(N), ...
## ======================================================================


## ----------------------------------------------------------------------
## Code 5.2 & 5.3: R code for analysis of 'Daphnia magna' data.
## ----------------------------------------------------------------------
library( nlreg )

# Fit transform both sides (TBS) model with constant variance and 
# produce diagnostic plot

daphnia.nl <- 
   nlreg( log(time) ~ log(b1+(b2-b1)/(1+(conc/b4)^(2*b3))), 
          data = daphnia, 
          start = list(b1 = 0.05, b2 = 50, b3 = 1.5, b4 = 0.2) ) 
summary( daphnia.nl )

nlreg.diag.plots( daphnia.nl, which = 2 )

# TBS model fit, profile plot and 95% confidence intervals for hormesis 
# parameter

daphnia.h.nl <- 
   nlreg( log(time) ~ log(b1+(b2-b1)/(1+(((conc-x0)/b4)^2)^b3)), 
          start = list(b1 = 0.05, b2 = 50, b3 = 3, b4 = 0.1, x0 = 0), 
          data = daphnia, hoa = TRUE )
coef( daphnia.h.nl )

daphnia.h.prof <- profile( daphnia.h.nl, offset = x0 )

plot( daphnia.h.prof )
summary( daphnia.h.prof )
## ======================================================================


## ----------------------------------------------------------------------
## Code 5.4 & 5.5: R code for analysis of RIA data.
## ----------------------------------------------------------------------
library( nlreg )

# Fit nonlinear model to RIA data

ria.nl <- nlreg( formula = count ~ b1+(b2-b1)/(1+(conc/b4)^b3), 
                 weights = ~ (b1+(b2-b1)/(1+(conc/b4)^b3))^g, 
                 data = ria, hoa = TRUE, 
                 start = c(b1=1.6, b2=20, b3=2, b4=330, g=2) )
summary( ria.nl )

# Profile plots and profile pair sketches

ria.prof <- profile( ria.nl )
par( mar = c(1.5, 1, 2, 1) )
contour( ria.prof, alpha = 0.05, cl1 = "black", cl2 = "black", 
         lwd2 = 2 )

# Compute maximum adjusted profile likelihood estimates of variance 
# parameters

ria.mpl <- mpl( ria.nl )
summary( ria.mpl )
## ======================================================================


## ----------------------------------------------------------------------
## Code 5.6: R code used for analysis of leukemia data.
## ----------------------------------------------------------------------
# Function to compute minus log likelihood for example U ("exU") of 
# Cox and Snell (1981), including some censoring if appropriate

library( AA )
exU <- leuk

# Minus log likelihood function for Weibull model

llweibreg <- function(logb0 = 4, b1 = -1, tau = 0) {
  y <- log( exU$t[exU$status==1] )  # observed
  mu1 <- logb0 + b1 * (exU$x[exU$status==1] - mean(exU$x))
  t1 <- sum( tau + exp((y-mu1) * exp(-tau)) - (y-mu1)*exp(-tau) )
  c <- log( exU$t[exU$status==0] )  # censored
  mu2 <- logb0 + b1 * (exU$x[exU$status==0] - mean(exU$x))
  t2 <- sum( exp((c-mu2) * exp(-tau)) )
  t1 + t2 
}

# Maximum likelihood fit to original data (no censoring)

exU$status <- 1

library( stats4 )
( full.fit <- mle(llweibreg) )

# This is the loop for profiling

b1.v <- seq(-2, 0, length = 200)
logj2 <- b0.b1 <- profloglik.b1 <- vector( "numeric", 
                                           length = length(b1.v) )
for( i in seq(along = b1.v) ) {
  fit <- mle( llweibreg, fixed = list(b1=b1.v[i]) )
  logj2[i] <- - log( det(vcov(fit)) )
  b0.b1[i] <- coef(fit)[1]
  profloglik.b1[i] <- as.numeric( logLik(fit) ) 
}
 
# Higher order calculations

profloglik.b1 <- profloglik.b1 - max(profloglik.b1)
lp1 <- predict( smooth.spline(b1.v, profloglik.b1), b1.v, 1 )
rho <- 0.5 * ( log( det(vcov(full.fit)) ) + logj2 )
q <- lp1$y * exp(rho)
r.p <- sign(q) * sqrt( -2 * profloglik.b1 )
rstar <- r.p + (1/r.p) * log(q/r.p)

# Confidence intervals by spline smooths avoiding middle for r*

rstarshort <- c( rstar[1:94], rstar[106:200] )
b1.vshort <- c( b1.v[1:94], b1.v[106:200] )
predict( smooth.spline(rstarshort, b1.vshort),
         c(2.33, 1.96, 1.64, -1.64, -1.96, -2.33), 0 )$y
## ======================================================================


## ----------------------------------------------------------------------
## Code 5.7: R code for analysis of PET film data.
## ----------------------------------------------------------------------
make.bits <- function(psi, lam, d) {  
#  auxiliary function to save space below
  xl <- log( d$volt - 5 + exp(lam[4]) )
  x0 <- 4.9
  xl0 <- log( x0 - 5 + exp(lam[4]) )
  kap <- exp( lam[2] + lam[3]*log(d$volt) )
  kap0 <- exp( lam[2] + lam[3]*log(x0) )
  beta0 <- log(psi) + 3*log(10) - lam[1]*xl0 - lgamma(1 + 1/kap0)
  gam <- exp( beta0 + lam[1]*xl )
  y <- d$time
  z <- (y/gam)^kap
  list( z = z, gam = gam, kap = kap, xl = xl, 
        inv = lam[1] * exp(lam[4]) / (d$volt - 5 + exp(lam[4])) ) }

nlogL <- function(psi, lam, d) {
#  negative log likelihood
  aux <- make.bits(psi, lam, d)
  sum( aux$z - d$cens * log(aux$kap * aux$z / d$time) ) }

library( boot )

nlogL.init <- function(th, d)  nlogL( th[1], th[-1], d )    
fit <- nlm( nlogL.init, c(1, 0, 0, 0, 0), d = hirose )
hirose.init <- fit$estimate   # initial values for maximisation

make.V <- function(th, d) {
#  make matrix containing V
  aux <- make.bits(th[1], th[-1], d)  
  V <- matrix(NA, nrow(d), 5)
  V[,c(1,2,5)] <- cbind(1, aux$xl, aux$inv)
  V[,3] <- - log(d$time / aux$gam )
  V[,4] <- - log(d$time / aux$gam) * log(d$volt)
  V * d$time }

phi <- function(th, V, d) {
#  use V to compute local parametrization
  aux <- make.bits(th[1], th[-1], d)
  L <- (d$cens * (aux$kap-1) - aux$z * aux$kap) / d$time
  apply(L*V, 2, sum) }

library( fraser.reid )
hirose.fr <- fraser.reid( NULL, nlogL = nlogL, phi = phi, make.V = make.V, 
                          th.init = jitter(hirose.init), d = hirose )
plot.fr(hirose.fr)
lik.ci(hirose.fr)
## ======================================================================


============================
CHAPTER 6: Some Case Studies
============================


## ----------------------------------------------------------------------
## Code 6.1: R code for approximate conditional inference for the grazing 
##           data.
## ----------------------------------------------------------------------
library( AA )
fit <- glm( State ~ Block + Species + Height, family = binomial, 
            data = treePref )   
anova( fit, test = "Chisq" )
summary(fit) 

library( cond )
fit.cond <- cond( fit, offset = HeightSmall, to = 7 ) 
summary(fit.cond)
## ======================================================================


## ----------------------------------------------------------------------
## Code 6.2: R code for analysis of metsulfuron-methyl data set M2.
## ----------------------------------------------------------------------
library( nlreg )

# Fit TBS model with constant variance

M2.nl <- nlreg( 
  formula = log(area) ~ log(b1+(b2-b1)/(1+(dose/b4)^(2*b3))), 
  data = M2, start = c(b1 = 130, b2 = 2400, b3 = 2, b4 = 0.05), 
  hoa = TRUE )
summary( M2.nl )

# Graphical inspection of regression diagnostics

nlreg.diag.plots( M2.nl )
# ======================================================================


## ----------------------------------------------------------------------
## Code 6.3: Model fit for the metsulfuron-methyl data set M4.
## ----------------------------------------------------------------------
library( nlreg )

# Fit TBS model with POX variance (MLEs & MMPLEs)

M4.pox.nl <- 
   nlreg( formula = log(area) ~ log(b2/(1+(dose/b4)^(2*b3))), 
          weights = ~ (1+dose^g)^2, data = M4, hoa = TRUE,
          start = c(b2 = 1500, b3 = 1, b4 = 0.05, g = 0.3) )

M4.pox.mpl <- mpl( M4.pox.nl )
summary( M4.pox.mpl )
 
# Use of diagnostic plots

M4.pox.diag <- nlreg.diag( M4.pox.nl )
plot( M4.pox.diag, which = 3 )
plot( M4.pox.diag, which = 9, iden = TRUE )
## ======================================================================


## ----------------------------------------------------------------------
## Code 6.4: R code for analysis of herbicide data set C3.
## ----------------------------------------------------------------------
library( nlreg )

# Fit TBS model with EIV variance function

C3.nl <- 
  nlreg( formula = log(area) ~ log(b1+(b2-b1)/(1+(dose/b4)^(2*b3))),
         weights = ~ 1 + k * dose^g * 
                     ( (b2-b1)*(2*b3)*dose^(2*b3-1)/
                       (1+(dose/b4)^(2*b3))^2/b4^(2*b3)/ 
                       (b1+(b2-b1)/(1+(dose/b4)^(2*b3))) )^2 ,
         start = c(b1=2.2, b2=1700, b3=2.8, b4=0.28, g=2.7, k=1),
         data = C3, hoa = TRUE, trace = TRUE,
         control = list(x.tol=1e-12, rel.tol=1e-12, step.min=1e-12))
summary( C3.nl )
## ======================================================================


=========================
CHAPTER 7: Further Topics
=========================


## ----------------------------------------------------------------------
## Code 7.1: R code for higher order inference for $\s^2_\alpha$ for the 
##           voltage regulator data.
## ----------------------------------------------------------------------
n <- c( 8, 4, 7, 7, 4, 7, 8, 6, 6, 7 )
ybar <- c( 16.3375,  16.1625,  16.05357, 15.96429, 16.16875, 
           16.22143, 15.97812, 16,       16.29583, 16.09286 )
volty <- list( ss = 8.7614, ss.df = 162, ssr = 9.4808, ssr.df = 54, 
               ybar = ybar, n = n )

nlogL <- function(psi, lam,  d) {
#  compute negative log likelihood
  sig2a <- exp(psi)
  sig2e <- exp(lam[1])
  sig2 <- exp(lam[2])
  mu <- lam[3]
  L <- d$ss/sig2 + log(sig2)*d$ss.df +  
       d$ssr/(sig2 + 4*sig2e) + log(sig2 + 4*sig2e)*d$ssr.df 
  L/2 - sum( dnorm(d$ybar, mu, sqrt( sig2a + sig2e/d$n + sig2/(4*d$n) ), 
                   log = T) ) }

make.V <- function(th, d) {
#  outputs n x dim(th) matrix with rows containing Vs
  sig2a <- exp(th[1])
  sig2e <- exp(th[2])
  sig2 <- exp(th[3])
  mu <- th[4]
  V <- matrix(NA, 12, 4)
  V[1,] <- c(0, 0, d$ss/sig2, 0)
  V[2,] <- c(0, 4, 1, 0) * d$ssr / (sig2 + 4*sig2e)
  S <- sqrt( sig2a + sig2e/d$n + sig2/(4 * d$n) )
  z <- (d$ybar - mu)/S
  V[3:12,-4] <- cbind( z/(2*S), z/(2*S)/d$n, z/(2*S)/(4*d$n) )
  V[3:12,4] <- 1
  V }

phi <- function(th, V, d) {
#  outputs vector of length dim(th)
  sig2a <- exp(th[1])
  sig2e <- exp(th[2])
  sig2 <- exp(th[3])
  mu <- th[4]
  S <- sig2a + sig2e/d$n + sig2/(4 * d$n)
  z <- (d$ybar - mu)/S
  L <- - c( 1/sig2, 1/(sig2 + 4*sig2e), 2*z )/2
  apply(V*L, 2, sum) }

volt.init <- c(0, 0, 0, 16)
L <- function(th, d)  nlogL( th[1], th[-1], d )
volt.fit <- nlm( L, volt.init, d = volty )  
# get initial values for likelihood fit

library( fraser.reid )
volt.fr <-  fraser.reid( NULL, nlogL, phi, make.V, volt.fit$estimate,
                         d = volty )
volt.ci <- lik.ci( volt.fr )
## ======================================================================


## ----------------------------------------------------------------------
## Code 7.2: R code for inference on $\beta$ based on return data, under
##           heteroscedastic error model.
## ----------------------------------------------------------------------
library( boot )

nlogL <- function(psi, lam, d)
{
  sig2 <- exp( lam[1] + lam[2]*d$market )
  mu <- lam[3] + psi*d$market
  - sum( dnorm(d$acme, mu, sqrt(sig2), log = T) )
}

make.V <- function(th, d)
{  
  sig2 <- exp( th[2] + th[3]*d$market )
  mu <- th[4] + th[1]*d$market
  z <- d$acme - mu
  cbind( 1, d$market, z, z*d$market )
}

phi <- function(th, V, d)
{ 
  sig2 <- exp( th[2] + th[3]*d$market )
  mu <- th[4] + th[1]*d$market
  L <- - (d$acme - mu)/sig2
  apply( V*L, 2, sum )
}

return.init <- rep(1, 4)
L <- function(th, d)  nlogL( th[1], th[-1], d )
return.fit <- nlm( L, return.init, d = acme )   
                                          # overall ML fit

library( fraser.reid )
return.fr1 <-  fraser.reid( NULL, nlogL, phi, make.V, 
                            return.fit$estimate, d = acme )
return.ci <- lik.ci( return.fr1 )
plot.fr( return.fr1 )
## ======================================================================


## ----------------------------------------------------------------------
## Code 7.3: R code for inference on air-conditioning data.  
##           The function "fitgamma" fits separate and common gamma 
##           models to a data frame of failure times and group number, 
##           and uses the "gamma.shape" function in the "MASS" package. 
##           The function "bootgamma" gives the parametric bootstrap 
##           distribution of the likelihood ratio statistic.
## ----------------------------------------------------------------------
fitgamma <- function(data) { 
  k <- nlevels( factor(data$number) ) 
  ydot <- tily <- n <- betahat <- muhat <- rep(0, k)
# Fit separate gamma models
  for(i in 1:k) {
    y <- data$time[data$number==i]
    ydot[i] <- sum(y)
    n[i] <- length(y)
    tily[i] <- sum(log(y))
    muhat[i] <- mean(y)
    betahat[i] <- as.numeric( 
                    gamma.shape( glm(y ~ 1, family = Gamma) )[1] ) 
  }
# Fit a common shape parameter	
    betall <- as.numeric( gamma.shape(
      glm(data$time ~ factor(data$number) - 1, family = Gamma) )[1] )
# Likelihood ratio statistic	
    loglikrat <- 2 * ( sum( - n*lgamma(betahat) + 
                            n*betahat*log(betahat/muhat) +
                            betahat*tily - (betahat/muhat)*ydot ) - 
                       sum( - n*lgamma(betall) + 
                            n*betall*log(betall/muhat) +
                            betall*tily - (betall/muhat)*ydot ) )
    list( k = k, n = n, ydot = ydot, tily = tily, muhat = muhat, 
          betahat = betahat, betall = betall, loglikrat = loglikrat ) 
}        

bootgamma <- function(data, fullfit, B) {  
  wboot <- rep(0, B)
  num <- data$number
  k <- fullfit$k
  n <- fullfit$n
  betall <- fullfit$betall
  muhat <- fullfit$muhat
  for(b in 1:B) {
    ystar <- as.list(rep(0, k))
   for(i in 1:k) {
      ystar[[i]] <- rgamma( n[i], shape = betall, 
                            scale = muhat[i]/betall )        
    }  
    aircondstar <- data.frame( time = unlist(ystar), number = num )
    wboot[b] <- fitgamma(data = aircondstar)$loglikrat 
  }
  wboot
}             
# ======================================================================


## ----------------------------------------------------------------------
## Code 7.4: Extraction of components of the "summary.cond" object for 
##           computing confidence intervals for the permutation 
##           distribution.
## ----------------------------------------------------------------------
library( AA )
boothb   # data frame with response and group indicator

booth.glm <- glm( Group ~ y, family = binomial, data = boothb )

library( cond )
booth.cond <- cond.glm( booth.glm, offset = y )
booth.summary <- summary( booth.cond, test = 0 )
# attributes( booth.summary ) lists the components of this object
booth.summary$signif.tests$stats[1:4, 2]
              
# Now we iterate "cond" over plausible range of lower limits

delta.low <- seq(1.5, 2, length = 30)  
conf.low <- matrix(0, nrow = length(delta.low), ncol = 3)
              
for( i in 1:length(delta.low) ) {
  shift <- boothb$y - delta.low[i]*boothb$Group
  shift.cond <- cond( glm( boothb$Group ~ shift, family = binomial ), 
                      offset = shift )
  shift.cond.summary <- summary( shift.cond, 
                                 test = 0 )$signif.tests$stats 
  conf.low[i,1] <- as.numeric( shift.cond.summary[4,2] )
  conf.low[i,2] <- as.numeric( shift.cond.summary[3,2] )
  conf.low[i,3] <- as.numeric( shift.cond.summary[1,2] ) 
}

predict( smooth.spline(conf.low[,1], delta.low), 0.025 )$y
predict( smooth.spline(conf.low[,2], delta.low), 0.025 )$y
predict( smooth.spline(conf.low[,3], delta.low), 0.025 )$y
## ======================================================================
