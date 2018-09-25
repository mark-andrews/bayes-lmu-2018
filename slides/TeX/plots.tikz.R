library(tikzDevice)
# TikZ device plot of Bernoulli likelihood.

n <- 250
m <- 139

beamer.parms = list(paperwidth   = 364.19536/72,  # converts `pt` to `in`
                    paperheight  = 273.14662/72,
                    textwidth    = 307.28987/72,
                    textheight   = 269.14662/72)

##########################################################################
tikz(file='binomial.likelihood.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

par(cex.lab=1.0,
    cex.axis=1.0,
    mar = c(4,2.5,0,0) )

y.label = "$\\theta^m  (1-\\theta)^{n-m}$"
x.label = "$\\theta$"

curve(x^m*(1-x)^(n-m), 
      from = 0, 
      to = 1, 
      n=1001,
      ylab='',
      xlab='',
      yaxt='n',
      axes=F)

x.ticks <- seq(0, 1, 0.1)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext(x.label, side=1, line=3)

max.limit <- 1.0 * ((m/n)^m*(1-(m/n))^(n-m))
y.ticks <- seq(0, max.limit, length.out = 5)
axis(2, at=y.ticks, labels=rep('',5), col.axis="black", las=1)
mtext(y.label, side=2, line=1.5)

dev.off()
###############################################################

###############################################################
plot.beta <- function(filename, alpha, beta, y.label='$\\mathrm{P}(\\theta)$'){
	tikz(file=filename,
	     standAlone=F,
	     width = beamer.parms$paperwidth * 0.9, 
	     height = beamer.parms$paperheight * 0.75)

	par(cex.lab=1.0,
	    cex.axis=1.0,
	    mar = c(4,2.5,0,0) )

	x.label = "$\\theta$"

	curve(dbeta(x, alpha, beta), 
	      ylab='',
	      xlab='',       
	      from=0, 
	      to=1, 
	      yaxt='n',
	      xaxt='n',
	      axes=F,
	      n=1001)

	x.ticks <- seq(0, 1, 0.1)
	axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
	mtext(x.label, side=1, line=3)

	max.limit <- dbeta(alpha/(alpha+beta), alpha, beta)
	y.ticks <- seq(0, max.limit, length.out = 5)
	axis(2, at=y.ticks, labels=rep('',5), col.axis="black", las=1)
	mtext(y.label, side=2, line=1.5)

	dev.off()
}
##################################################################

source('../code/utils/beta.hpd.R')

beta.plot.hpd <- function(filename, alpha, beta, y.label='$\\mathrm{P}(\\theta)$'){
  
  tikz(file=filename,
       standAlone=F,
       width = beamer.parms$paperwidth * 0.9, 
       height = beamer.parms$paperheight * 0.75)
  
  par(cex.lab=1.0,
      cex.axis=1.0,
      mar = c(4,2.5,0,0) )
  
  HPD <- beta.hpd.interval(alpha, beta)
  
  hpd.interval <- HPD$hpd.interval
  p.star <- HPD$p.star
  
  curve(dbeta(x, alpha, beta), 
        ylab='',
        xlab='',
        from=0, 
        to=1,
        yaxt='n',
        xaxt='n',
        axes=F,
        n=1001)
  
  y.label = "$\\mathrm{P}(\\theta)$"
  x.label = "$\\theta$"
  
  segments(hpd.interval[1], p.star, hpd.interval[2], p.star, lwd=1)
  
  x.ticks <- seq(0, 1, 0.1)
  axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
  mtext(x.label, side=1, line=3)
  
  max.limit <- dbeta((alpha-1)/(alpha+beta-2), alpha, beta)
  y.ticks <- seq(0, max.limit, length.out = 5)
  axis(2, at=y.ticks, labels=rep('',5), col.axis="black", las=1)
  mtext(y.label, side=2, line=1.5)
  
  dev.off()
}
##################################################################
plot.beta('beta.3.5.distribution.tex', 3, 5)
plot.beta('beta.140.112.distribution.tex', 140, 112, '$\\mathrm{P}(\\theta\\vert D)$')

beta.plot.hpd('beta.3.5.hpd.tex', 3, 5, '$\\mathrm{P}(\\theta\\vert D)$')
beta.plot.hpd('beta.140.112.hpd.tex', 140, 112, '$\\mathrm{P}(\\theta\\vert D)$')
##################################################################

##################################################################
tikz(file='gamma.distributions.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.75, 
     height = beamer.parms$paperheight * 0.6)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,0,0) )

colors <- c('chocolate',
            'coral4',
            'cadetblue4')

curve(dgamma(x, shape=2.0, scale=10),
      from = 0.0,
      to = 100,
      col=colors[1],
      ylim=c(0, 0.04),
      xlab='',
      ylab='',
      axes=F,
      n=1001)

curve(dgamma(x, shape=3.0, scale=12),
      from = 0.0,
      to = 100,
      col=colors[2],
      n=1001, add=T)

curve(dgamma(x, shape=1.0, scale=100),
      from = 0.0,
      to = 100,
      col=colors[3],
      n=1001, add=T)

legend(50, 0.035, 
       legend=c('$\\kappa=2$, $\\theta=10$', '$\\kappa=3$, $\\theta=12$', '$\\kappa=1.0$, $\\theta=100$'),
       lty=1, 
       col=colors, 
       bty='n')

x.ticks <- seq(0, 100, 10)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=3)

y.ticks <- seq(0, 0.04, by = 0.01)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(\\lambda)$', side=2, line=4.0)

dev.off()
########################################################################

tikz(file='poisson.likelihood.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.75, 
     height = beamer.parms$paperheight * 0.65)

par(cex.lab=1.0,
    cex.axis=1.0,
    mar = c(3.0,2.5,0,0) )

n <- 10
S <- 107
curve((exp(-n*x)*x^S)/(exp(-n*(S/n))*(S/n)^S),
      from = 5.0,
      to = 20,
      xlab='',
      ylab='',
      col='salmon4',
      axes=F,
      n=1001)

x.ticks <- seq(5, 20, 5)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=2.0)

y.ticks <- seq(0, 1.0, by = 0.2)
axis(2, at=y.ticks, labels=rep('', length(y.ticks)), col.axis="black", las=1)
mtext('$\\textrm{P}(D\\vert\\lambda)$', side=2, line=1.5)

dev.off()


##################################################################
tikz(file='gamma.posterior.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(3,4,1,0) )

curve(dgamma(x, shape=108, scale=1/(10+1/100)),
      from = 0.0,
      to = 20,
      col='salmon4',
      #ylim=c(0, 0.04),
      xlab='',
      ylab='',
      axes=F,
      n=1001)

x.ticks <- seq(0, 100, 10)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=2)

y.ticks <- seq(0, 0.4, by = 0.1)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(\\lambda\\vert D)$', side=2, line=3.0)

dev.off()

source('../code/utils/gamma.R')

tikz(file='gamma.posterior.hpd.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(3,4,1,0) )

kappa <- 108
theta <- 1/(10+1/100)

HPD <- gamma.hpd.interval(kappa, theta)
hpd.interval <- HPD$hpd.interval
p.star <- HPD$p.star

curve(dgamma(x, shape=kappa, scale=theta),
      from = 0.0,
      to = 20,
      col='salmon4',
      #ylim=c(0, 0.04),
      xlab='',
      ylab='',
      axes=F,
      n=1001)

segments(hpd.interval[1], p.star, hpd.interval[2], p.star, lwd=1, col='salmon4')

x.ticks <- seq(0, 100, 10)
axis(1, at=x.ticks, labels=x.ticks, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=2)

y.ticks <- seq(0, 0.4, by = 0.1)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(\\lambda\\vert D)$', side=2, line=3.0)

dev.off()

# Also shows the hpd interval.
gamma.plot.hpd.2 <- function(alpha, beta){
  
  kappa <- 108
  theta <- 1/(10+1/100)
  
  HPD <- gamma.hpd.interval(kappa, theta)
  hpd.interval <- HPD$hpd.interval
  p.star <- HPD$p.star
  
  x <- seq(hpd.interval[1],
           hpd.interval[2],
           length.out = 1001)
  
  y <- dgamma(x, kappa, scale=theta)
  x <- c(hpd.interval[1], x, hpd.interval[2])
  y <- c(0, y, 0)
  
  plot.new()
  plot(NULL, NULL,
       ylab=expression("P" * (lambda)),
       xlab=expression(lambda),
       xlim=c(0,20),
       ylim=c(0, 1.1 * dgamma((kappa-1)*theta, kappa, scale=theta))
  )
  
  polygon(x, y, col="lavender", lty=0)
  
  curve(dgamma(x, kappa, scale=theta), 
        ylab='',
        xlab='',
        from=0, 
        to=20, 
        n=1001,
        add=T)
  
  #segments(hpd.interval[1], p.star, hpd.interval[2], p.star, lwd=1)
  
}

tikz(file='poisson.posterior.predictive.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

V <- posterior.predictive(x)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,1,0) )

x.new <- seq(0, 25)
p <- dnbinom(x.new, size=V$r, prob = 1-V$q)
plot(x.new, 
     p, 
     ylim=c(0, 0.11),
     ylab='',
     xlab='',
     type='h', 
     col='chocolate',
     axes=F)

axis(1, at=x.new, labels=x.new, col.axis="black", las=1)
mtext('$k$', side=1, line=3)

y.ticks <- seq(0, 0.1, by=0.02)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(k\\vert D)$', side=2, line=4)

dev.off()



############################################################

tikz(file='poisson.distribution.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.75, 
     height = beamer.parms$paperheight * 0.6)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,1,0) )

x.new <- seq(0, 15)
p <- dpois(x.new, lambda=5)
plot(x.new, 
     p, 
     ylim=c(0, 0.18),
     ylab='',
     xlab='',
     type='h', 
     col='chocolate',
     axes=F)

axis(1, at=x.new, labels=x.new, col.axis="black", las=1)
mtext('$\\lambda$', side=1, line=3)

y.ticks <- seq(0, 0.15, by=0.05)
axis(2, at=y.ticks, labels=y.ticks, col.axis="black", las=1)
mtext('$\\textrm{P}(x=k\\vert \\lambda)$', side=2, line=4)

dev.off()

source('../code/utils/normal.data.R')
source('../code/utils/normal.R')

tikz(file='normal.likelihood.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.75, 
     height = beamer.parms$paperheight * 0.6)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,1,0) )

contour(mu.lim, 
        sigma.lim, 
        xlim=c(3.2, 5),
        ylim=c(0.8, 2.1),
        L,
        nlevels = 5,
        axes=F,
        ylab='',
        xlab='')

x.new <- seq(3.2, 5, length.out=7)
axis(1, at=x.new, labels=round(x.new, 2), col.axis="black", las=1)
mtext('$\\theta$', side=1, line=3)

y.ticks <- seq(0.8, 2.1, length.out=5)
axis(2, at=y.ticks, labels=round(y.ticks, 2), col.axis="black", las=1)
mtext('$\\sigma$', side=2, line=4)

dev.off()

########################################################################
source('../code/utils/normal.R')
source('../code/utils/invgamma.R')
tikz(file='scaled.inv.chisq.tex',
     standAlone=F,
     width = beamer.parms$paperwidth * 0.9, 
     height = beamer.parms$paperheight * 0.75)

colors <- c('chocolate',
            'coral4',
            'salmon',
            'cadetblue4')
n <- 1e3
nu <- c(5, 10, 2, 1)
tau <- c(1.1, 2, 2, 1)

par(cex.lab=1.2,
    cex.axis=1.0,
    mar = c(4,5,1,0) )

curve(dinvchisq(x, nu[1], tau[1]), from=0.0, to=10.0, col=colors[1], axes=F, n=n, xlab='', ylab='')
curve(dinvchisq(x, nu[2], tau[2]), from=0.0, to=10.0, n=n, col=colors[2], add=T)
curve(dinvchisq(x, nu[3], tau[3]), from=0.0, to=10.0, col=colors[3], add=T, n=n)
curve(dinvchisq(x, nu[4], tau[4]), from=0.0, to=10.0, col=colors[4], add=T, n=n)

labels <- c()
for (k in c(1,2,3,4)){
  labels <- c(labels, paste("$\\nu = ", nu[k], "$,  $\\tau = ", tau[k], "$", sep=''))
}
legend(5, 0.5, legend=labels, col=colors, lwd=1, border=NA, bty='n')


x.new <- seq(0, 10, length.out=10)
axis(1, at=x.new, labels=round(x.new, 2), col.axis="black", las=1)
mtext('$x$', side=1, line=3)

y.ticks <- seq(0.0, 0.6, by=0.1)
axis(2, at=y.ticks, labels=round(y.ticks, 2), col.axis="black", las=1)
mtext('$P(x\\vert \\nu, \\tau^2)$', side=2, line=4)

dev.off()
