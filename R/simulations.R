# Simulation analysis
# Starting on 2 dimensions with prior on covariance matrix: 
# 1) IW(v,L), v=d+1=3 always. Options for L : 
#   - Identity, Diagonal with sample var, Pre-scale data and use Identity
# 2) SIW... regular one ? 
# 3) SS with correlations from IW and variaces with IG prior
# 
# Scenarios: 
#   dimension=2, s1 = .01, 1, 100 and s2 = .01, 1, 100, rho=0, .5, .99
# Output: 
#   density plots, 
# how is the right way to "scale back" posterior inference on the pre-scale data case?

library(plyr)
library(reshape2)
library(ggplot2)
library(mnormt)
library(rstan)
set_cppo(mode = "fast")

# set up the parallel for plyr functions
# parallel <- require(doMC, quietly=TRUE)
# if(parallel){
#   registerDoMC(8)
# }

# Compile the model objects 
source('R/modelcode.R')
m_iw  <- stan_model(model_code=sim.iw)
m_siw <- stan_model(model_code=sim.siw)
m_ss  <- stan_model(model_code=sim.ssig)
#save(m_iw, m_siw, m_ssig, file='data/models_cpp.Rdata')

# functions to run stan model
runstan.sim <- function(d, it = 1500,ch=3, w=500, prm=NULL,prior=NULL) {
  if (prior == 'iw')  mod<- m_iw
  if (prior == 'iw2')  mod<- m_iw
  if (prior == 'siw') mod<- m_siw
  if (prior == 'ss.ig')  mod<- m_ss
  K <- ncol(d[,-c(1:5)])
  
  if (prior == 'iw2') R <- diag(apply(d[,-c(1:5)],2,var))
  if (prior != 'iw2') R <- diag(rep(1,K))

  dat = list(y = d[,-c(1:5)], N = nrow(d), R = R, k=K, mu0 = rep(0,K))
  sampling(object=mod, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
}
printresult <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}
getiter <- function(xx) {
  attributes(xx)$stan_args[[2]]$iter
}

simula <- function(size, data,...) {
  prms <- c('s1', 's2', 'rho')
  simdata <- subset(data, ns==size)                       
  
#  mod_iw <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='iw',...)  
  mod_iw2 <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='iw2',...)
  #  mod_iw3 <-  dlply(simdata, .(sim,r,s,ns),runstan.sim, prm=prms, prior='iw3',...)
  mod_siw <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='siw',...)                        
  mod_ssig <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='ss.ig',...)                        
#data.frame(prior='iw',ldply(mod_iw, printresult)),
  
res.df <- rbind(   data.frame(prior='iw2',ldply(mod_siw, printresult)),
                   data.frame(prior='siw',ldply(mod_siw, printresult)),
                   data.frame(prior='ss.ig',ldply(mod_ssig, printresult)) )
  out <- list(res.df, mod_iw2, mod_siw, mod_ssig)
  names(out) <- c('res', 'iw2', 'siw', 'ss.ig')
  return(out)
}

load('data/simdata.Rdata')
res_size10d2 <- simula(size=10, data=simdata.2,it=600,w=100)
save(res_size10d2, file='data/sims_n10_d2.Rdata')

# res_size50d2 <- simula(size=50, data=data2)
# save(res_size50d2, file='../data/sims_n50_d2.Rdata')
# res_size250d2 <- simula(size=250, data=data2)
# save(res_size250d2, file='../data/sims_n250_d2.Rdata')

# ============
# For testing...
d <- subset(simdata.2, s1==.01 & s2==.01 & r==.99 & ns==10 & sim==1)
dat = list(y = d[,-c(1:5)], N = nrow(d), R = diag(2), k=2, mu0 = rep(0,2))
toy <- sampling(object=m_siw, data = dat)
, pars=c('s1','s2','rho'), iter = 1100, chains = 3, warmup=100)
toy
# =========

