# RESULTS

library(ggplot2)
library(plyr)
library(reshape2)
library(rstan)
library(xtable)

load('data/sims_n10_d2.Rdata')
res <- res_size10d2$res

# Gelman diag
summary(res$Rhat)

diag <- subset(res, Rhat > 2 & param!='lp__', c(1:7,14,17))
with(diag, table(prior, param))
with(subset(diag, prior=='ss.ig'), table(s1,s2))


f1 <- function(x) {
  # function to extract samples form the stanfit object
  xx <- extract(x,permuted=F)
  xx <- dcast(melt(xx), iterations+chains~parameters)
  colnames(xx)[4:6] <- paste('sim',c('rho','s1','s2'),sep='.')
  return(xx)
}


samples <- ldply(res_size10d2[-1], function(x) ldply(x,f1) )



pdf('simplot.pdf')
d <- subset(res, param=='s2')
d$r <- factor(d$r, levels=c(0,0.5,0.99), labels=paste('rho',c(0,0.5,0.99),sep='=='))
#d$s <- factor(d$s, levels=c(0.01,0.1,1,10,100), labels=paste('sigma',c(0.01,0.1,1,10,100),sep='=='))
qplot(data=d ,x=s2, y=mean,color=factor(s1),xlab='True SD2', ylab='Posterior Mean SD2') + 
  facet_grid(facets=r~prior,scales='free',labeller=label_parsed) + geom_abline(1) + theme(legend.position= 'bottom') + scale_x_log10() + scale_y_log10()

d <- subset(res, param=='rho')
d$rx <- d$r + runif(nrow(d),-.1,.1)
#d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$s1 <- factor(d$s1, levels=c(0.01,0.1,1,100), labels=paste('sigma',c(0.01,0.1,1,100),sep='=='))
d$s2 <- factor(d$s2, levels=c(0.01,0.1,1,100), labels=paste('sigma',c(0.01,0.1,1,100),sep='=='))
qplot(data=d ,x=rx, y=mean,color=prior,shape=prior,xlab='True Correlation', ylab='Correlation Posterior Mean') + 
  facet_grid(facets=s2~s1,scales='free',labeller=label_parsed) + geom_abline(1) + theme(legend.position= 'bottom') + scale_x_continuous(breaks=c(0,0.5,1))

d <- subset(samples, r==0)
qplot(data=d,x=sim.s2,shape=factor(sim),geom='density',color=.id, main='posterior density for s2 when rho=0') +
  facet_grid(facets=s1~s2,scales='free') + theme(legend.position= 'bottom') + scale_x_log10()

d <- subset(samples, r==0.5)
qplot(data=d,x=sim.s2,shape=factor(sim),geom='density',color=.id, main='posterior density for s2 when rho=0.5') +
  facet_grid(facets=s1~s2,scales='free') + theme(legend.position= 'bottom') + scale_x_log10()

d <- subset(samples, r==0.99)
qplot(data=d,x=sim.s2,shape=factor(sim),geom='density',color=.id, main='posterior density for s2 when rho=0.99') +
  facet_grid(facets=s1~s2,scales='free') + theme(legend.position= 'bottom') + scale_x_log10()

sd.s2 <- ddply(samples,.(.id,r,s1,s2,sim), summarise, sd2 =log(sd(sim.s2)) )
d <- dcast(sd.s2, sim+r+s1+s2 ~ .id)
qplot(data=d, x=iw2,y=siw, main='are we cheating?')+facet_wrap(s1~s2,scales='free') + geom_abline(1) 
dev.off()

