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
m_iw   <- stan_model(model_code=sim.iw)
m_sciw <- stan_model(model_code=sim.sciw)
m_siw  <- stan_model(model_code=sim.siw)
m_ss.igiw   <- stan_model(model_code=sim.ssigiw)
m_ss.iglkj   <- stan_model(model_code=sim.ssiglkj)
m_ss.lnlkj   <- stan_model(model_code=sim.sslnlkj)
m_ss.lniw   <- stan_model(model_code=sim.sslniw)

#save(m_iw, m_siw, m_ssig, file='data/models_cpp.Rdata')

# functions to run stan model
runstan.sim <- function(d, it = 1500,ch=3, w=500, prm=NULL,prior=NULL) {
# Priors: 
  # iw   :inverse wishart with identity and df=dim+1
  # iw2  :inverse wishart with sample variance on diag and df=dim+1
  # iw3  :inverse wishart with prescaled data
  # siw  : Scaled inverse wishart, with lognormal for scaled parameters

  # ss.igiw  : Separation strategy using correlations from IW and variance as IG
  # ss.iglkj : Separation strategy using correlations from LKJ and variance as IG
  # ss.lniw  : Separation strategy using correlations from IW and variance as LN
  # ss.lnlkj : Separation strategy using correlations from LKJ and variance as LN
  
  if (prior == 'iw')   mod<- m_iw
  if (prior == 'iw2')  mod<- m_iw
  if (prior == 'iw3')  mod<- m_sciw
  if (prior == 'siw')  mod<- m_siw
  if (prior == 'ss.igiw')  mod<- m_ss.igiw
  if (prior == 'ss.iglkj')  mod<- m_ss.iglkj
  if (prior == 'ss.lniw')  mod<- m_ss.lniw
  if (prior == 'ss.lnlkj')  mod<- m_ss.lnlkj

  K <- ncol(d[,-c(1:5)])
  
  if (prior == 'iw3') {
    aux <- scale(d[,-c(1:5)], center=FALSE)
    sample.sd <- attributes(aux)[[3]]
    ds <- cbind(d[,1:5], aux)  
    dat = list(y = ds[,-c(1:5)], N = nrow(ds), R = diag(2), k=2, mu0 = rep(0,2), samplesd=sample.sd)
  }
  if (prior == 'iw2') R <- diag(apply(d[,-c(1:5)],2,var))
  
  if (prior != 'iw2') R <- diag(rep(1,K))
  if (prior != 'iw3') dat = list(y = d[,-c(1:5)], N = nrow(d), R = R, k=K, mu0 = rep(0,K))
  
  out <- sampling(object=mod, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
  x <- printresult(out) 
  gd <- max(x$Rhat)
  it2 <- it + 1000; w2 <- w+500; j<-1
  stay.in.loop <- gd > 1.1
  
  while (stay.in.loop) {
    out <- sampling(object=mod, data = dat,pars=prm, iter = it2, chains = ch, warmup=w2)
    x  <- printresult(out)
    gd <- max(x$Rhat, na.rm=T)
    j <- j + 1
    it2 <- it2 + 1500
    w2 <- 1500
    stay.in.loop <- (gd > 1.1) & (j < 6)
  }
  return(out)
}

printresult <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}

simula <- function(size, data,priorlist=c('iw2', 'siw', 'ss.ig', 'iw3'), ...) {
  prms <- c('s1', 's2', 'rho')
  simdata <- subset(data, ns==size)                       
  out <- llply(priorlist, function(x) dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior=x,...)  )
  names(out) <- priorlist
  out$res <- ldply(out, function(x) ldply(x, printresult) )
  return(out)
}

#simula.old <- function(size, data,...) {
  prms <- c('s1', 's2', 'rho')
  simdata <- subset(data, ns==size)                       
  
  #  mod_iw <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='iw',...)  
  mod_iw2 <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='iw2',...)
  mod_siw <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='siw',...)                        
  mod_ssig <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='ss.ig',...)                        
  mod_iw3 <-  dlply(simdata, .(sim,r,s1,s2,ns),runstan.sim, prm=prms, prior='iw3',...)
  #data.frame(prior='iw',ldply(mod_iw, printresult)),
  
  res.df <- rbind(   data.frame(prior='iw2',ldply(mod_siw, printresult)),
                     data.frame(prior='siw',ldply(mod_siw, printresult)),
                     data.frame(prior='ss.ig',ldply(mod_ssig, printresult)) )
  out <- list(res.df, mod_iw2, mod_siw, mod_ssig, mod_iw3)
  names(out) <- c('res', 'iw2', 'siw', 'ss.ig', 'iw3')
  return(out)
}

load('data/simdata.Rdata')

# simulation for SS, 4 priors total
simSS_n10d2 <- simula(size=10, data=simdata.2, priorlist=c('siw','ss.igiw','ss.iglkj','ss.lniw', 'ss.lnlkj'),it=2000,w=500)
save(simSS_n10d2, file='data/simSS_n10d2.Rdata')

# simulations for IW variants
res_size10d2 <- simula(size=10, data=simdata.2, priorlist=c('siw', 'iw', 'iw2', 'iw3'),it=600,w=100)
save(res_size10d2, file='data/sims_n10_d2.Rdata')

# res_size50d2 <- simula(size=50, data=data2)
# save(res_size50d2, file='../data/sims_n50_d2.Rdata')
# res_size250d2 <- simula(size=250, data=data2)
# save(res_size250d2, file='../data/sims_n250_d2.Rdata')


#######################################################################
#######################################################################
#######################################################################
# ============
# For testing...
d <- subset(simdata.2, s1==.01 & s2==.01 & r==0 & ns==10 & sim==1)
aux <- scale(d[,-c(1:5)], center=FALSE)
sample.sd <- attributes(aux)[[3]]
ds <- cbind(d[,1:5], aux)
dat = list(y = ds[,-c(1:5)], N = nrow(ds), R = diag(2), k=2, mu0 = rep(0,2), samplesd=sample.sd)
toy <- sampling(object=m_siw, data = dat)
, pars=c('s1','s2','rho'), iter = 1100, chains = 3, warmup=100)
toy
x <- extract(toy, permuted=FALSE)
xx <- dcast(melt(x), iterations+chains~parameters)
xx$uns1 <- sqrt(xx[, 'Sigma[1,1]']  )
qplot(data=melt(xx[,c('s1','uns1')]),x=value,,geom='density',color=variable) + scale_x_log10()
toy2 <- simula(size=10, data=d, priorlist=c('ss.iglkj','ss.igiw'),it=2000,w=500)

# =========

