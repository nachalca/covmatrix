
# Covariance Matrix project: focus on IW
# Create data for simulations. 

# Scenario: Multivariate normal data, with 0 mean and different cov matrix 
# dimension=2 
# standard dev: s1 = .01, 1, 100 and s2 = .01, 1, 100 
# correlations: rho=0, .5, .99

library(plyr)
library(mnormt)

# We first simulate with Variance = 1 and then re-scale to obtain other data sets. 
simdat <- function(ns,r,dim) {
  # r: is the common rho value for all dimensions
  sigma <- diag(dim); l <- lower.tri(sigma); u <- upper.tri(sigma)
  sigma[l] <- r; sigma[u]<- r
  mu <- rep(0,dim)
  data.frame(rmnorm(ns, mean=mu, varcov=sigma))
}  

# Simulate 5 data sets for each combination of r,size keeping separate the dimension
prm <- expand.grid(r= c(0,.5,.99), ns=c(10,50,250))

#set.seed(1234)

s50 <-  rdply(5,mdply(prm, simdat, dim=50))
colnames(s50)[1] <- 'sim'
s2 <- s50[,1:5]
#s10 <- s100[,1:13]

# rescale each data set, ss represent the new standar deviation
# there are 450 data set in total
m <- matrix(1:16, ncol=4)
ss <- expand.grid(s1=c(.01,.1,1,100), s2=c(.01,.1,1,100))
ss <- ss[lower.tri(m, diag=TRUE),]

rescale    <- function(s1,s2, dx) data.frame(dx[,1:3],X1=dx$X1*s1,X2=dx$X2*s2)
simdata.2    <- mdply(ss, rescale, dx=s2)

# simdata.10   <- mdply(ss, rescale, dx=s10)
# simdata.100  <- mdply(ss, rescale, dx=s100)
# number of data set: 
# xx <- ddply(simdata.2,.(ns,s1,s2,r,sim), summarise, m = sd(X1))
# xx <- adply(simdata.2[,1:5], .margins=1, paste, colapse='.')
# prod(laply(xx, dim)[,1])

# Save data 
save(simdata.2, file='data/simdata.Rdata')
