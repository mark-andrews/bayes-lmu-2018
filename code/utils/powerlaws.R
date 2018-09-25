library(VGAM)

likelihood.f <- function(lambda, K, xm){
  if (lambda > xm){
    y <- lambda^-K
  } else {
    y <- 0 
  }
  return(y)
}

likelihood.f <- Vectorize(likelihood.f)

K <- 5
xm <- 20

par(mfrow=c(2,1))

# Likelihood function
curve(likelihood.f(x, K=K, xm=xm),
      xlim=c(0, 100),
      xlab=expression(lambda),
      ylab=expression(P*(x*'|'*lambda)))
      
curve(dpareto(x, scale=xm, shape=K-1),
      xlim=c(0, 100),
      xlab=expression(lambda),
      ylab=expression(P*(lambda*'|'*x)))