

# intro to Bayes workshop functions


# this function is not intuitively named but is the optional stopping demo 
# you can call it with just t.updated() or change the parameters


t.updated <- function(initial.n=10, increment=1, plot=TRUE, like=TRUE){
  iterations <- 1000
  data <- rnorm(initial.n)
  if (increment==1) index <- 0:iterations+initial.n
  else index = 1:1001
  p.vector <- vector('numeric', iterations)
  lr.vector <- vector('numeric', iterations)
  p.initial <- t.test(data)$p.val
  mu.initial <- 0.5 # t.test(data)$est[[1]]
  lr.initial <- LR.t(mu1=t.test(data)$est[[1]], mu2=0, t.test(data)$est[[1]],
                     SE=t.test(data)$est[[1]]/t.test(data)$stat[[1]], df= length(data)-1, independent=FALSE)[[3]]
  for (i in 1:iterations) {
    data <- c(data, rnorm(increment))
    p.vector[i] <- t.test(data)$p.val
    lr.vector[i] <- LR.t(mu1=mu.initial, mu2=0, t.test(data)$est[[1]],
                         SE=t.test(data)$est[[1]]/t.test(data)$stat[[1]], df= length(data)-1, independent=FALSE)[[3]]
  }
  p.vect <- as.vector(c(p.initial, p.vector))
  lr.vect <- as.vector(c(lr.initial, lr.vector))
  plot(index, p.vect, ylim=c(0,1), ylab='p value', xlab='Optional stopping', pch=20)
  abline(h=.05, col='red')
  if(like ==TRUE) {
    par(new=TRUE, mar=c(4,4,4,4))
    plot(index, log(lr.vect), ylab="", xlab='Optional stopping', col='green', pch=20, axes=FALSE)
    abline(h=1, col='blue')
    axis(side=4)
    mtext("Log(LR) in favour of null",side=4,col="black", line=2.5)
  }
  output <- list(p.vect,lr.vect)
  output
}






# functions from the Baguley (2012) Serious Stats book
# other functions from the book
# source('http://www2.ntupsychology.net/seriousstats/SeriousStatsAllfunctions.txt')


# updated 29 03 15

# This file contains all functions from Thom Baguley's (2012) Serious Stats book and some additional functions in a single file.


# functions from Chapter 6 of Serious Stats

ssc.R <- function(R, N, q = 1) {
  # Ezekiel (1929) small sample correction for r
  R.out <- (1 - (((1 - R^2) * (N - 1))/(N - q - 1)))^0.5
  if (R < 0) 
    R.out <- R.out * -1
  R.out
}


rz.ci <- function(r, N, conf.level = 0.95) {
  zr.se <- 1/(N - 3)^0.5
  moe <- qnorm(1 - (1 - conf.level)/2) * zr.se
  zu <- atanh(r) + moe
  zl <- atanh(r) - moe
  tanh(c(zl, zu))
} 


r.ind.ci <- function(r1, r2, n1, n2=n1, conf.level = 0.95) {
  L1 <- rz.ci(r1, n1, conf.level = conf.level)[1]
  U1 <- rz.ci(r1, n1, conf.level = conf.level)[2]
  L2 <- rz.ci(r2, n2, conf.level = conf.level)[1]
  U2 <- rz.ci(r2, n2, conf.level = conf.level)[2]
  lower <- r1 - r2 - ((r1 - L1)^2 + (U2 - r2)^2)^0.5
  upper <- r1 - r2 + ((U1 - r1)^2 + (r2 - L2)^2)^0.5
  c(lower, upper)
} 


# functions from Chapter 7 of Serious Stats

smd.d <- function(m1, m2, s1, s2, n1, n2 = n1) {
  num <- m2 - m1
  denom <- (((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2))^0.5
  num/denom
}

smd.g <- function(m1, m2, s1, s2, n1, n2 = n1) {
  num <- m2 - m1
  denom <- (((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))^0.5
  num/denom
}

smd.unb <- function(m1, m2, s1, s2, n1, n2 = n1) {
  num <- m2 - m1
  denom <- (((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))^0.5
  num/denom * (n1 + n2 - 3)/(n1 + n2 - 2.25)
} 

t.to.d.rm <- function(t, n) {
  t * (1/n)^0.5 * (n/(n - 1))^0.5
}

t.to.g.rm <- function(t, n) {
  t * (1/n)^0.5
}

t.to.g.unb.rm <- function(t, n) {
  g <- t.to.g.rm(t, n)
  g * (n - 3)/(n - 2.25)
} 

t.to.d.rm <- function(t, n) {
  t * (1/n)^0.5 * (n/(n - 1))^0.5
}

t.to.g.rm <- function(t, n) {
  t * (1/n)^0.5
}

t.to.g.unb.rm <- function(t, n) {
  g <- t.to.g.rm(t, n)
  g * (n - 3)/(n - 2.25)
} 

g.to.r <- function(g, n1, n2 = n1) {
  N <- n1 + n2
  g/(g^2 + (N^2 - 2 * N)/(n1 * n2))^0.5
}

r.to.g <- function(r, n1, n2 = n1) {
  N <- n1 + n2
  r/((1 - r^2) * ((n1 * n2)/(N^2 - 2 * N)))^0.5
}

r.true <- function(r.obs, rxx) {
  r.obs/sqrt(rxx * r.yy)
} 


# functions from Chapter 10 of Serious Stats

boot.rlm <- function(data, indices, maxit=20) {
  data <- data[indices,]
  rlm(y ~ x, method='MM', data=data, maxit=maxit)$coef[[2]]
}



# functions from Chapter 11 of Serious Stats

binom.lik <- function(successes, trials, accuracy = 1e-05, plot = FALSE) {
  P <- seq(0, 1, accuracy)
  like <- P^successes * (1 - P)^(trials - successes)
  max.like <- max(like)
  like <- like/max.like
  P.hat <- P[which(like == max(like))]
  b8.i <- min(which(like >= 1/8))
  e8.i <- max(which(like >= 1/8))
  b8.i <- min(which(like >= 1/8))
  e8.i <- max(which(like >= 1/8))
  b32.i <- min(which(like >= 1/32))
  e32.i <- max(which(like >= 1/32))
  b32.i <- min(which(like >= 1/32))
  e32.i <- max(which(like >= 1/32))
  b8 <- P[b8.i]
  e8 <- P[e8.i]
  b32 <- P[b32.i]
  e32 <- P[e32.i]
  output <- list(`P-hat` = P.hat, `1/8 likelihood interval` = c(b8, e8), 
                 `1/8 lower` = b8, `1/8 upper` = e8, `1/32 likelihood interval` = c(b32, 
                                                                                    e32), `1/32 lower` = b32, `1/32 upper` = e32)
  if (plot == TRUE) {
    curve((x^successes * (1 - x)^(trials - successes)/max.like), xlim = c(0, 
                                                                          1), xlab = expression(italic("P")), ylab = "Likelihood")
    segments(b8, 1/8, e8, 1/8, lwd = 0.2)
    segments(b32, 1/32, e32, 1/32, lwd = 0.2)
  }
  output
} 


pois.lik <- function(count, plot = FALSE) {
  range <- count * 10 + 10
  lambda <- seq(0, range, range/10000)
  like <- lambda^count * exp(-1 * lambda)
  max.like <- max(like)
  like <- like/max.like
  pois.mle <- count
  b8.i <- min(which(like >= 1/8))
  e8.i <- max(which(like >= 1/8))
  b8.i <- min(which(like >= 1/8))
  e8.i <- max(which(like >= 1/8))
  b8 <- lambda[b8.i]
  e8 <- lambda[e8.i]
  b32.i <- min(which(like >= 1/32))
  e32.i <- max(which(like >= 1/32))
  b32.i <- min(which(like >= 1/32))
  e32.i <- max(which(like >= 1/32))
  b32 <- lambda[b32.i]
  e32 <- lambda[e32.i]
  output <- list(`1/8 likelihood interval` = c(b8, e8), `1/8 lower` = b8, 
                 `1/8 upper` = e8, `1/32 likelihood interval` = c(b32, e32), `1/32 lower` = b32, 
                 `1/32 upper` = e32)
  if (plot == TRUE) {
    curve(x^count * exp(-1 * x)/max.like, xlim = c(0, count * 2 + 7), xlab = expression(lambda), 
          ylab = "Likelihood")
    segments(b8, 1/8, e8, 1/8, lwd = 0.2)
    segments(b32, 1/32, e32, 1/32, lwd = 0.2)
  }
  output
} 


t.lik <- function(t, df, independent = TRUE, adjust = FALSE) {
  if (independent == TRUE) 
    N = df + 2
  else N = df + 1
  if (adjust == TRUE) 
    N = df
  like <- (1 + t^2/df)^-((N)/2)
  like
}

LR.t <- function(mu1, mu2, mu.obs, SE, df, independent = TRUE, adjust = FALSE) {
  lik1 <- t.lik((mu1 - mu.obs)/SE, df, independent = independent, adjust = adjust)
  lik2 <- t.lik((mu2 - mu.obs)/SE, df, independent = independent, adjust = adjust)
  LR.mu1 <- lik1/lik2
  LR.mu2 <- 1/LR.mu1
  return <- list("Likelihood ratios", LR.mu1 = LR.mu1, LR.mu2 = LR.mu2)
  return
} 

t.lik.int <- function(mean, SE, df, independent = TRUE, adjust = FALSE, 
                      accuracy = 1e-04, plot = FALSE) {
  if (independent == TRUE) 
    N = df + 2
  else N = df + 1
  SDest <- SE * N^0.5
  if (adjust == TRUE) 
    N <- df
  mu <- seq(mean - SDest, mean + SDest, accuracy)
  like <- (1 + ((mean - mu)/SE)^2/df)^-((N)/2)
  max.like <- max(like)
  M.hat <- mean
  b8.i <- min(which(like >= 1/8))
  e8.i <- max(which(like >= 1/8))
  b8 <- mu[b8.i]
  e8 <- mu[e8.i]
  b32.i <- min(which(like >= 1/32))
  e32.i <- max(which(like >= 1/32))
  b32 <- mu[b32.i]
  e32 <- mu[e32.i]
  if (plot == TRUE) {
    curve((1 + ((mean - x)/SE)^2/df)^-(N/2), xlim = c(-SDest + mean, SDest + 
                                                        mean), xlab = expression(mu), ylab = "Likelihood")
    segments(b8, 1/8, e8, 1/8, lwd = 0.2)
    segments(b32, 1/32, e32, 1/32, lwd = 0.2)
  }
  output <- list(`M-hat` = M.hat, `1/8 likelihood interval` = c(b8, e8), 
                 `1/8 lower` = b8, `1/8 upper` = e8, `1/32 likelihood interval` = c(b32, 
                                                                                    e32), `1/32 lower` = b32, `1/32 upper` = e32)
  output
} 


Bayes.norm.1s <- function(M.obs, SD.obs, mu.prior, sigma.prior, 
                          N, probability = 0.95, ssc = TRUE, plot = FALSE) {
  if (ssc == TRUE) 
    SD.obs <- SD.obs * (1 + 20/N^2)
  SE.adj <- SD.obs/N^0.5
  precision.prior <- 1/sigma.prior^2
  precision.lik <- 1/SE.adj^2
  precision.post <- precision.prior + precision.lik
  mean.post <- (precision.prior/precision.post) * mu.prior + (precision.lik/precision.post) * 
    M.obs
  sigma.post <- (1/precision.post)^0.5
  lower <- mean.post + qnorm((1 - probability)/2) * sigma.post
  upper <- mean.post - qnorm((1 - probability)/2) * sigma.post
  post.prob.int <- c(lower, upper)
  min.m <- min(M.obs, mu.prior, mean.post)
  max.m <- max(M.obs, mu.prior, mean.post)
  max.s <- max(SE.adj, sigma.prior, sigma.post)
  y.max <- max(dnorm(mu.prior, mu.prior, sigma.prior), dnorm(M.obs, M.obs, 
                                                             SE.adj), dnorm(mean.post, mean.post, sigma.post))
  if (plot == TRUE) {
    curve(dnorm(x, mean.post, sigma.post), lty = 3, xlim = c(min.m - 3.75 * 
                                                               max.s, max.m + 3.75 * max.s), yaxt = "n", ylab = NA, xlab = expression(mu))
    curve(dnorm(x, mu.prior, sigma.prior), lwd = 0.3, add = TRUE)
    curve(dnorm(x, M.obs, SE.adj), add = TRUE)
    legend(max.m + max.s, y.max, legend = c("Prior", "Likelihood", "Posterior"), 
           lty = c(1, 1, 3), lwd = c(0.3, 1, 1))
  }
  output <- list(`Posterior mean` = mean.post, `Posterior SD` = sigma.post, 
                 `Probability interval` = post.prob.int, Probability = probability, lower = lower, 
                 upper = upper)
  output
} 

Bayes.norm.2s <- function(M.diff, SD.1, SD.2, diff.prior, sigma.prior, 
                          n1, n2, probability = 0.95, ssc = TRUE, plot = FALSE) {
  if (ssc == TRUE) {
    SD.1 <- SD.1 * (1 + 20/n1^2)
    SD.2 <- SD.2 * (1 + 20/n2^2)
  }
  SS.pooled <- SD.1^2 * (n1 - 1) + SD.2^2 * (n2 - 1)
  SE.adj <- (SS.pooled/(n1 + n2 - 2) * (1/n1 + 1/n2))^0.5
  precision.prior <- 1/sigma.prior^2
  precision.lik <- 1/SE.adj^2
  precision.post <- precision.prior + precision.lik
  mean.post <- (precision.prior/precision.post) * diff.prior + (precision.lik/precision.post) * 
    M.diff
  sigma.post <- (1/precision.post)^0.5
  lower <- mean.post + qnorm((1 - probability)/2) * sigma.post
  upper <- mean.post - qnorm((1 - probability)/2) * sigma.post
  post.prob.int <- c(lower, upper)
  min.m <- min(M.diff, diff.prior, mean.post)
  max.m <- max(M.diff, diff.prior, mean.post)
  max.s <- max(SE.adj, sigma.prior, sigma.post)
  y.max <- max(dnorm(diff.prior, diff.prior, sigma.prior), dnorm(M.diff, M.diff, 
                                                                 SE.adj), dnorm(mean.post, mean.post, sigma.post))
  if (plot == TRUE) {
    curve(dnorm(x, mean.post, sigma.post), lty = 3, xlim = c(min.m - 3.75 * 
                                                               max.s, max.m + 3.75 * max.s), yaxt = "n", ylab = NA, xlab = expression(mu))
    curve(dnorm(x, diff.prior, sigma.prior), lwd = 0.3, add = TRUE)
    curve(dnorm(x, M.diff, SE.adj), add = TRUE)
    legend(max.m + max.s, y.max, legend = c("Prior", "Likelihood", "Posterior"), 
           lty = c(1, 1, 3), lwd = c(0.3, 1, 1))
  }
  output <- list(`Posterior mean` = mean.post, `Posterior SD` = sigma.post, 
                 `Probability interval` = post.prob.int, Probability = probability, lower = lower, 
                 upper = upper)
  output
} 


unit.prior.Bf.1s <- function(t, N, scale.factor = 1) {
  df <- N - 1
  numerator <- (1 + t^2/df)^-((df + 1)/2)
  denom.1 <- (1 + N * scale.factor^2)^-0.5
  denom.2 <- (1 + t^2/(df * (1 + N * scale.factor^2)))^-((df + 1)/2)
  Bf <- numerator/(denom.1 * denom.2)
  output <- list(`Bayes factor for H0` = Bf, `Bayes factor for H1` = 1/Bf)
  output
}

unit.prior.Bf.2s <- function(t, n1, n2 = n1, scale.factor = 1) {
  df <- n1 + n2 - 2
  N <- (1/n1 + 1/n2)^-1
  numerator <- (1 + t^2/df)^-((df + 1)/2)
  denom.1 <- (1 + N * scale.factor^2)^-0.5
  denom.2 <- (1 + t^2/(df * (1 + N * scale.factor^2)))^-((df + 1)/2)
  Bf <- numerator/(denom.1 * denom.2)
  output <- list(`Bayes factor for H0` = Bf, `Bayes factor for H1` = 1/Bf)
  output
}

# Note the following functions written by Danny Kaye and Thom Baguley

# the following functions written by Danny Kaye and Thom Baguley
# adding in scale factor r and changing defaults for r
# note that the default scale factor in Serious Stats (2012) was r = 1 as per Rouder et al. (2009)
# the new functions change the defaults to match those of the new BayesFactor package (Morey & Rouder, 2013)

denOneSample<-function(g, N, t, r) {
  e <- exp(1)
  nu <- N - 1
  ret <- (1+N*g*r^2)^(-1/2)*(1+t^2/((1+N*g*r^2)*nu))^(-(nu+1)/2)*(2*pi)^(-1/2)*g^(-3/2)*e^(-1/(2*g))
  ret
}

denTwoSample<-function(g, N, t, r, nu) {
  e <- exp(1)
  ret <- (1+N*g*r^2)^(-1/2)*(1+t^2/((1+N*g*r^2)*nu))^(-(nu+1)/2)*(2*pi)^(-1/2)*g^(-3/2)*e^(-1/(2*g))
  ret
}

numTwoSample<-function(N, t, nu) {
  ret <- (1+(t^2)/nu)^(-(nu+1)/2)
  ret
}

numOneSample<-function(N, t) {
  nu <- N-1
  ret <- (1+(t^2)/nu)^(-(nu+1)/2)
  ret
}

crit<-function(N, t, c) {
  ret <- c+1
  while(ret > c)
  {
    t <- t + 0.01
    ret <- BfOneSample(N, t)
  }
  t
}

JZS.prior.Bf.1s <- function(t, N, r = 0.5){
  a <- integrate(denOneSample,0,Inf,N,t,r)
  Bf <- numOneSample(N, t)/ a$value
  output <- list('Bayes factor for H0' = Bf, 'Bayes factor for H1' = 1/Bf, 'Scale factor' = r)
  output
}

JZS.prior.Bf.2s <- function(t, n1, n2 = n1, r = sqrt(2)/2) {
  N <- n1 * n2 / (n1 + n2)
  nu = n1 + n2 - 2
  a <- integrate(denTwoSample,0,Inf,N,t,r,nu)
  Bf <- numTwoSample(N, t, nu)/ a$value
  output <- list('Bayes factor for H0' = Bf, 'Bayes factor for H1' = 1/Bf, 'Scale factor' = r)
  output
}


# functions from Chapter 12 of Serious Stats

r.partial <- function(r.ab, r.ac, r.bc) {
  numerator <- r.ab - r.ac * r.bc
  denominator <- ((1 - r.ac^2) * (1 - r.bc^2))^0.5
  output <- numerator/denominator
  output
}



# functions from Chapter 14 of Serious Stats

# Note that the following function are prediction equations for specific models and not generally useful

p.main <- function(present, future){3.366-0.2727*present-0.0294*future}

p.int <- function(present, future){5.7935-1.2309*present-0.1586*future+0.0506*present*future}

# This is a generic function for interaction/profile plots of 2x2 designs from cell means

plot.2by2 <- function(A1B1, A1B2, A2B1, A2B2, group.names, legend = TRUE, 
                      leg.loc = NULL, factor.labels = c("Factor A", "Factor B"), swap = FALSE, 
                      ylab = NULL, main = NULL) {
  group.means <- c(A1B1, A2B1, A1B2, A2B2)
  if (missing(ylab)) 
    ylab <- expression(italic(DV))
  if (swap == TRUE) {
    group.names <- list(group.names[[2]], group.names[[1]])
    group.means <- c(A1B1, A1B2, A2B1, A2B2)
    factor.labels <- c(factor.labels[2], factor.labels[1])
  }
  plot(group.means, pch = NA, ylim = c(min(group.means) * 0.95, max(group.means) * 
                                         1.025), xlim = c(0.8, 2.2), ylab = ylab, xaxt = "n", xlab = factor.labels[1], 
       main = main)
  points(group.means[1:2], pch = 21)
  points(group.means[3:4], pch = 19)
  axis(side = 1, at = c(1:2), labels = group.names[[1]])
  lines(group.means[1:2], lwd = 0.6, lty = 2)
  lines(group.means[3:4], lwd = 0.6)
  if (missing(leg.loc)) 
    leg.loc <- c(1, max(group.means))
  if (legend == TRUE) 
    legend(leg.loc[1], leg.loc[2], legend = group.names[[2]], title = factor.labels[2], 
           lty = c(3, 1))
} 



# functions from Chapter 15 of Serious Stats

nu.prime <- function(means, sds, ns, weights = 1) {
  sum(weights^2 * sds^2/ns)^2/sum((weights^4 * sds^4)/(ns^2 * (ns - 1)))
}


wts.diff <- function(weights.A, weights.B, rescale = FALSE) {
  sd.A <- (sum(weights.A^2)/length(weights.A))^0.5
  sd.B <- (sum(weights.B^2)/length(weights.B))^0.5
  contrast.weights <- weights.A/sd.A - weights.B/sd.B
  if (rescale == TRUE) 
    contrast.weights <- 2/sum(abs(contrast.weights)) * contrast.weights
  contrast.weights
}



# functions from Chapter 16 of Serious Stats

# the plot functions are taken from Baguley (2012) within-subject CI paper
# the latest versions are added at the end of this file


# functions from Chapter 18 of Serious Stats

pb.sim <- function(model) {
  sim <- simulate(model)
  fixef(refit(model,sim))
}


# function from Online Supplements for Serious Stats

mi.to.mitools <- function(imputed.data.from.mi, m = imputed.dat@m) {
  # depends on mi and mitools
  data.list <- as.list(1:m)
  for (i in 1:m) data.list[[i]] <- mi.data.frame(imputed.data.from.mi, m = i)
  mitools.list <- imputationList(data.list)
  mitools.list
} 


# functions from http://seriousstats.wordpress.com/2012/02/05/comparing-correlations/

rho.rxy.rxz <- function(rxy, rxz, ryz) {
  num <- (ryz - 1/2*rxy*rxz) * (1 - rxy^2 - rxz^2 - ryz^2) + ryz^3
  den <- (1 - rxy^2) * (1 - rxz^2)
  num/den
}


r.dol.ci <- function(r12, r13, r23, n, conf.level = 0.95) {
  L1 <- rz.ci(r12, n, conf.level = conf.level)[1]
  U1 <- rz.ci(r12, n, conf.level = conf.level)[2]
  L2 <- rz.ci(r13, n, conf.level = conf.level)[1]
  U2 <- rz.ci(r13, n, conf.level = conf.level)[2]
  rho.r12.r13 <- rho.rxy.rxz(r12, r13, r23)
  lower <- r12 - r13 - ((r12 - L1)^2 + (U2 - r13)^2 - 2 * rho.r12.r13 * (r12 - L1) * (U2 - r13))^0.5
  upper <- r12 - r13 + ((U1 - r12)^2 + (r13 - L2)^2 - 2 * rho.r12.r13 * (U1 - r12) * (r13 - L2))^0.5
  c(lower, upper)
} 


rho.rab.rcd <- function(rab, rac, rad, rbc, rbd, rcd) {
  num <- 1/2 * rab*rcd * (rac^2 + rad^2 + rbc^2 + rbd^2) + rac*rbd + rad*rbc - (rab*rac*rad + rab*rbc*rbd + rac*rbc*rcd + rad*rbd*rcd)
  den <- (1 - rab^2) * (1 - rcd^2)
  num/den
}

r.dnol.ci <- function(r12, r13, r14, r23, r24, r34, n12, n34=n12, conf.level = 0.95) {
  L1 <- rz.ci(r12, n12, conf.level = conf.level)[1]
  U1 <- rz.ci(r12, n12, conf.level = conf.level)[2]
  L2 <- rz.ci(r34, n34, conf.level = conf.level)[1]
  U2 <- rz.ci(r34, n34, conf.level = conf.level)[2]
  rho.r12.r34 <- rho.rab.rcd(r12, r13, r14, r23, r24, r34)
  lower <- r12 - r34 - ((r12 - L1)^2 + (U2 - r34)^2 - 2 * rho.r12.r34 * (r12 - L1) * (U2 - r34))^0.5
  upper <- r12 - r34 + ((U1 - r12)^2 + (r34 - L2)^2 - 2 * rho.r12.r34 * (U1 - r12) * (r34 - L2))^0.5
  c(lower, upper)
} 


# functions from Baguley (2012) Calculating and graphing within-subject confidence intervals for ANOVA
# see http://psychologicalstatistics.blogspot.com/2011/10/calculating-and-graphing-within-subject.html
# they have now been updated slightly - see URL <to be added>


# note: intend to rewrite function to remove dependence on nlme
lm.ci <- function(data.frame, conf.level = 0.95, difference = FALSE) {
  #loftus-masson within-subject CIs
  k = ncol(data.frame)
  n <- nrow(data.frame)
  df.stack <- stack(data.frame)
  require(nlme)
  parts <- rep(1:n, k)
  root.ms.error <- lme(values ~ 0 + ind, random = ~1 | parts, cbind(parts, 
                                                                    df.stack))[[6]]
  detach(package:nlme)
  mean.mat <- matrix(, k, 1)
  ci.mat <- matrix(, k, 2)
  if (difference == TRUE) 
    diff.factor = 2^0.5/2
  else diff.factor = 1
  moe <- root.ms.error/n^0.5 * qt(1 - (1 - conf.level)/2, (n - 1) * (k - 
                                                                       1)) * diff.factor
  for (i in 1:k) mean.mat[i, ] <- colMeans(data.frame[i])
  for (i in 1:k) {
    ci.mat[i, 1] <- mean.mat[i] - moe
    ci.mat[i, 2] <- mean.mat[i] + moe
  }
  dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
  ci.mat
}

bs.ci <- function(data.frame, conf.level = 0.95, difference = FALSE) {
  # between-subject CIs
  k = ncol(data.frame)
  n <- nrow(data.frame)
  df.stack <- stack(data.frame)
  group.means <- colMeans(data.frame, na.rm = TRUE)
  if (difference == TRUE) 
    ci.mat <- (confint(lm(values ~ 0 + ind, df.stack)) - group.means) * 
    2^0.5/2 + group.means
  else ci.mat <- confint(lm(values ~ 0 + ind, df.stack))
  dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
  ci.mat
}

cm.ci.u <- function(data.frame, conf.level = 0.95) {
  #cousineau uncorrected within-subject CIs
  k = ncol(data.frame)
  n <- nrow(data.frame)
  df.stack <- stack(data.frame)
  index <- rep(1:n, k)
  p.means <- tapply(df.stack$values, index, mean)
  norm.df <- data.frame - p.means + (sum(data.frame)/(n * k))
  output <- matrix(, k, 2)
  dimnames(output) <- list(names(data.frame), c("lower", "upper"))
  for (i in 1:k) output[i, ] <- t.test(norm.df[i], conf.level = conf.level)$conf.int[1:2]
  output
}

cm.ci <- function(data.frame, conf.level = 0.95, difference = TRUE) {
  #cousineau-morey within-subject CIs
  k = ncol(data.frame)
  if (difference == TRUE) 
    diff.factor = 2^0.5/2
  else diff.factor = 1
  n <- nrow(data.frame)
  df.stack <- stack(data.frame)
  index <- rep(1:n, k)
  p.means <- tapply(df.stack$values, index, mean)
  norm.df <- data.frame - p.means + (sum(data.frame)/(n * k))
  t.mat <- matrix(, k, 1)
  mean.mat <- matrix(, k, 1)
  for (i in 1:k) t.mat[i, ] <- t.test(norm.df[i])$statistic[1]
  for (i in 1:k) mean.mat[i, ] <- colMeans(norm.df[i])
  c.factor <- (k/(k - 1))^0.5
  moe.mat <- mean.mat/t.mat * qt(1 - (1 - conf.level)/2, n - 1) * c.factor * 
    diff.factor
  ci.mat <- matrix(, k, 2)
  dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
  for (i in 1:k) {
    ci.mat[i, 1] <- mean.mat[i] - moe.mat[i]
    ci.mat[i, 2] <- mean.mat[i] + moe.mat[i]
  }
  ci.mat
}

ml.ci <- function(data.frame, conf.level = 0.95, cov.matrix = "unstructured") {
  # CI based on multilevel model with covariance matrix unstructured or
  #   constrained to compound symmetry
  k = ncol(data.frame)
  n.parts <- nrow(data.frame)
  data.long <- reshape(data.frame, idvar = "id", direction = "long", varying = 1:k, 
                       v.names = "dv")
  require(nlme)
  if (cov.matrix == "comp.symm") 
    ml.mod <- lme(dv ~ 0 + factor(time), random = ~1 | id, na.action = na.omit, 
                  data.long)
  else ml.mod <- lme(dv ~ 0 + factor(time), random = ~0 + factor(time) | 
                       id, na.action = na.omit, data.long)
  detach(package:nlme)
  require(gmodels)
  ci.mat <- ci(ml.mod, confidence = conf.level)[, 2:3]
  detach(package:gmodels)
  ci.mat
}

plot.wsci <- function(data.frame, conf.level = 0.95, type = "cm", 
                      difference = TRUE, cov.matrix = "unstructured", xlab = NULL, ylab = NULL, 
                      level.labels = NULL, main = NULL, pch.cex = 1.2, text.cex = 1.2, pch = 21,
                      ylim = c(min.y, max.y), line.width= c(1.5, 0), grid = FALSE) {
  # plot within-subject CIs by various methods
  k = ncol(data.frame)
  if (type == "cm") 
    ci.mat <- cm.ci(data.frame, conf.level, difference = difference)
  if (type == "uncorrected") 
    ci.mat <- cm.ci.u(data.frame, conf.level)
  if (type == "lm") 
    ci.mat <- lm.ci(data.frame, conf.level, difference = difference)
  if (type == "bs") 
    ci.mat <- bs.ci(data.frame, conf.level, difference = difference)
  if (type == "ml") 
    ci.mat <- ml.ci(data.frame, conf.level, cov.matrix = cov.matrix)
  moe.y <- max(ci.mat) - min(ci.mat)
  min.y <- min(ci.mat) - moe.y/3
  max.y <- max(ci.mat) + moe.y/3
  means <- colMeans(data.frame, na.rm = TRUE)
  if (missing(xlab)) 
    xlab <- "levels"
  if (missing(ylab)) 
    ylab <- "Confidence interval for mean"
  if (missing(level.labels) == FALSE) 
    names(data.frame) <- level.labels
  plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, k + 0.3), xlab = xlab, 
       ylab = ylab, main = main, cex.lab = text.cex)
  if (grid == TRUE) 
    grid()
  points(means, pch = pch, bg = "black", cex = pch.cex)
  index <- 1:k
  segments(index, ci.mat[, 1], index, ci.mat[, 2], lwd = line.width[1])
  segments(index - 0.02, ci.mat[, 1], index + 0.02, ci.mat[, 1], lwd = line.width[2])
  segments(index - 0.02, ci.mat[, 2], index + 0.02, ci.mat[, 2], lwd = line.width[2])
  axis(1, 1:k, labels = names(data.frame))
}

two.tiered.ci <- function(data.frame, conf.level = 0.95, cov.matrix = "unstructured", 
                          difference = TRUE, level.labels = NULL, xlab = NULL, ylab = NULL, main = NULL, 
                          pch = 19, pch.cex = 1.4, text.cex = 1.2, ylim = c(min.y, max.y), grid = FALSE,
                          line.width = c(1.5, 1.5), tier.width = 0) {
  # plot two tiered CI with ml approach for outer tier and cm approach for inner tier
  k = ncol(data.frame)
  ci.outer <- ml.ci(data.frame, conf.level = conf.level, cov.matrix = cov.matrix)
  moe.y <- max(ci.outer) - min(ci.outer)
  min.y <- min(ci.outer) - moe.y/3
  max.y <- max(ci.outer) + moe.y/3
  means <- colMeans(data.frame)
  if (missing(xlab)) 
    xlab <- "levels"
  if (missing(ylab)) 
    ylab <- "Confidence interval for mean"
  if (missing(level.labels) == FALSE) 
    names(data.frame) <- level.labels
  plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, k + 0.3), xlab = xlab, 
       ylab = ylab, main = main, cex.lab = text.cex)
  if (grid == TRUE) 
    grid()
  points(means, pch = pch, bg = "black", cex = pch.cex)
  index <- 1:k
  segments(index, ci.outer[, 1], index, ci.outer[, 2], lwd = line.width[1])
  axis(1, index, labels = names(data.frame))
  row.names(ci.outer) <- names(data.frame)
  if (cov.matrix == "comp.symm") 
    ci.inner <- lm.ci(data.frame, conf.level, difference = difference)
  else
    ci.inner <- cm.ci(data.frame, conf.level, difference = difference)
  if(tier.width == 0) {
    segments(index - 0.02, ci.inner[, 1], index + 0.02, ci.inner[, 1], lwd = line.width[2])
    segments(index - 0.02, ci.inner[, 2], index + 0.02, ci.inner[, 2], lwd = line.width[2])
  }
  else
    segments(index, ci.inner[, 1], index, ci.inner[, 2], lwd = line.width[1]*(1+abs(tier.width)))
}


# Note groups must be coded 1 to n and not start at zero
cm.ci.mixed <- function(data.frame, group.var = "last", conf.level = 0.95, 
                        difference = TRUE) {
  #cousineau-morey within-subject CIs for mixed design
  k = ncol(data.frame)
  if (difference == TRUE) 
    diff.factor = 2^0.5/2
  else diff.factor = 1
  if (group.var == "last") {
    within <- 1:(k - 1)
    data.frame[k] <- unclass(data.frame[[k]])[1:nrow(data.frame)]
  }
  else {
    within <- 2:k
    data.frame[1] <- unclass(data.frame[[1]])[1:nrow(data.frame)]
  }
  if (group.var == "last") 
    n.groups <- nlevels(factor(data.frame[[k]]))
  else n.groups <- nlevels(factor(data.frame[[1]]))
  c.factor <- ((k - 1)/(k - 2))^0.5
  ci.list <- vector("list", n.groups)
  for (i in 1:n.groups) {
    if (group.var == "last") 
      data <- subset(data.frame, data.frame[k] == i)[within]
    else data <- subset(data.frame, data.frame[1] == i)[within]
    p.means <- colMeans(as.data.frame(t(data)))
    norm.dat <- data - p.means + mean(p.means)
    t.mat <- matrix(, k - 1, 1)
    mean.mat <- matrix(, k - 1, 1)
    ci.mat <- matrix(, k - 1, 2)
    for (j in 1:(k - 1)) t.mat[j, ] <- t.test(norm.dat[j])$statistic[1]
    for (j in 1:(k - 1)) mean.mat[j, ] <- colMeans(norm.dat[j])
    n <- nrow(data)
    moe.mat <- mean.mat/t.mat * qt(1 - (1 - conf.level)/2, n - 1) * c.factor * 
      diff.factor
    for (j in 1:(k - 1)) {
      ci.mat[j, 1] <- mean.mat[j] - moe.mat[j]
      ci.mat[j, 2] <- mean.mat[j] + moe.mat[j]
    }
    dimnames(ci.mat) <- list(names(data), c("lower", "upper"))
    ci.list[[i]] <- ci.mat
  }
  ci.list
}

ml.ci.mixed <- function(data.frame, group.var = "last", conf.level = 0.95, 
                        cov.matrix = "unstructured") {
  # mixed CI based on multilevel model unstructured or constrained to compound symmetry
  k = ncol(data.frame)
  n.parts <- nrow(data.frame)
  if (group.var == "last") 
    within <- 1:(k - 1)
  else within <- 2:k
  data.long <- reshape(data.frame, idvar = "id", direction = "long", varying = within, 
                       v.names = "dv")
  group <- factor(data.long[[1]])
  n.groups <- nlevels(group)
  require(nlme)
  if (cov.matrix == "within.group.cs") 
    ml.mod <- lme(dv ~ 0 + factor(time):factor(group), random = ~ 0 + 
                    factor(time) | id, na.action = na.omit, cbind(group, data.long))
  if (cov.matrix == "unstructured") 
    ml.mod <- lme(dv ~ 0 + factor(time):factor(group), random = list(id = pdDiag(form = ~0 + 
                                                                                   factor(group)), id = ~0 + factor(time)), na.action = na.omit, 
                  cbind(group, data.long))
  if (cov.matrix == "comp.symm") 
    ml.mod <- lme(dv ~ 0 + factor(time):factor(group), random = ~1 | 
                    id, na.action = na.omit, cbind(group, data.long))
  require(gmodels)
  detach(package:nlme)
  ci.mat <- ci(ml.mod, confidence = conf.level)[, 2:3]
  group.id <- rep(1:n.groups, rep(k - 1, n.groups))
  output <- cbind(group.id, ci.mat)
  detach(package:gmodels)
  output
}


# note: intend to rewrite function to add tier.width argument
two.tiered.mixed <- function(data.frame, group.var = "last", conf.level = 0.95, 
                             cov.matrix = "unstructured", difference = TRUE, lines = FALSE, level.labels = NULL, 
                             xlab = NULL, ylab = NULL, main = NULL, pch = c(21:25, 1:3), pch.cex = 1.2, 
                             pch.col = c(3:9), text.cex = 1.2, ylim = c(min.y, max.y), grid = FALSE, jitter = NULL,
                             group.labels = c(1:8), leg.loc = c(1, min(ci.outer)), line.width= c(1.25, 1.5, 1),
                             tier.width = 0) {
  # plot two tiered mixed CI with ml approach for outer tier and cm
  #   approach for inner tier
  k = ncol(data.frame)
  if (group.var == "last") 
    n.groups <- nlevels(factor(data.frame[[k]]))
  else n.groups <- nlevels(factor(data.frame[[1]]))
  ci.outer <- ml.ci.mixed(data.frame, group.var = group.var, conf.level = conf.level, 
                          cov.matrix = cov.matrix)
  ci.inner <- cm.ci.mixed(data.frame, group.var = group.var, conf.level, 
                          difference = difference)
  ci.group.outer <- matrix(, k - 1, 2)
  ci.group.inner <- matrix(, k - 1, 2)
  index <- 1:(k - 1)
  moe.y <- max(ci.outer) - min(ci.outer)
  min.y <- min(ci.outer) - moe.y/3
  max.y <- max(ci.outer) + moe.y/3
  if (missing(xlab)) 
    xlab <- "levels"
  if (missing(ylab)) 
    ylab <- "Confidence interval for mean"
  if (missing(level.labels) == FALSE) 
    names(data.frame) <- level.labels
  plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, k - 1 + 0.3), xlab = xlab, 
       ylab = ylab, main = main, cex.lab = text.cex)
  if (grid == TRUE) 
    grid()
  axis(1, index, labels = rownames(ci.inner[[1]]))
  if (missing(jitter)) 
    jitter <- scale(1:n.groups, scale = FALSE)/(n.groups * 1.5)
  for (i in 1:n.groups) {
    ci.group.outer <- subset(as.data.frame(ci.outer), group.id == i)[2:3]
    ci.group.inner <- ci.inner[[i]]
    means <- (ci.group.inner[, 1] + ci.group.inner[, 2])/2
    segments(index + jitter[i], ci.group.outer[, 1], index + jitter[i], 
             ci.group.outer[, 2], lwd = line.width[1])
    if(tier.width == 0) {
      segments(index - 0.02 + jitter[i], ci.group.inner[, 1], index + 0.02 + 
                 jitter[i], ci.group.inner[, 1], lwd = line.width[2])
      segments(index - 0.02 + jitter[i], ci.group.inner[, 2], index + 0.02 + 
                 jitter[i], ci.group.inner[, 2], lwd = line.width[2])
    }
    else segments(index + jitter[i], ci.group.inner[, 1], index + jitter[i], 
                  ci.group.inner[, 2], lwd = line.width[1] * (1 + abs(tier.width)), col=pch.col[i])
    points(index + jitter[i], means, pch = pch[i], bg = pch.col[i], cex = pch.cex)
    if (lines == TRUE) 
      lines(index + jitter[i], means, lty = i + 1, lwd = line.width[3])
  }
  legend(leg.loc[1], leg.loc[2], legend = group.labels[1:n.groups], pch = pch[1:n.groups], 
         lty = 2:(n.groups + 1), pt.bg = pch.col[1:n.groups], bty = "n", horiz = TRUE, lwd = line.width[3])
} 


# these functions add between-subject versions of the functions
# see http://seriousstats.wordpress.com/2012/03/18/cis-for-anova/

# note input must number the factor levels 1 to J in order required for graphing
# also note input is in long format (standard for between-subjects data)

bsci <- function(data.frame, group.var=1, dv.var=2,  difference=TRUE, var.equal=FALSE, conf.level = 0.95) {
  data <- subset(data.frame, select=c(group.var, dv.var))
  fact <- factor(data[[1]])
  dv <- data[[2]]
  J <- nlevels(fact)
  N <- length(dv)
  ci.mat <- matrix(,J,3, dimnames=list(levels(fact), c('lower', 'mean', 'upper')))
  ci.mat[,2] <- tapply(dv, fact, mean)
  if(difference==TRUE) diff.factor= 2^.5/2 else diff.factor=1
  if(var.equal==TRUE) {
    moe <- summary(lm(dv ~ 0 + fact))$sigma/(N/J)^.5 * qt(1-(1-conf.level)/2,N-J) * diff.factor
    for(i in 1:J) {
      ci.mat[i,1] <- ci.mat[i,2] - moe
      ci.mat[i,3] <- ci.mat[i,2] + moe
    }
  }
  if(var.equal==FALSE) {
    for(i in 1:J) {
      group.dat <- subset(data, data[1]==levels(fact)[i])[[2]]
      n <- length(group.dat)
      moe <- sd(group.dat)/sqrt(n) * qt(1-(1-conf.level)/2,n-1) * diff.factor
      ci.mat[i,1] <- ci.mat[i,2] - moe
      ci.mat[i,3] <- ci.mat[i,2] + moe
    }
  }
  ci.mat
}



plot.bsci <- function(data.frame, group.var=1, dv.var=2,  difference=TRUE, var.equal=FALSE, conf.level = 0.95, xlab = NULL, ylab = NULL, level.labels = NULL, main = NULL, pch = 21, ylim = c(min.y, max.y), line.width= c(1.5, 0),grid=TRUE) {
  data <- subset(data.frame, select=c(group.var, dv.var))
  fact <- factor(data[[1]])
  dv <- data[[2]]
  J <- nlevels(fact)
  ci.mat <- bsci(data.frame=data.frame , group.var=group.var, dv.var=dv.var,  difference=difference, var.equal=var.equal, conf.level =conf.level)
  moe.y <- max(ci.mat) - min(ci.mat)
  min.y <- min(ci.mat) - moe.y/3
  max.y <- max(ci.mat) + moe.y/3
  if (missing(xlab)) 
    xlab <- "Groups"
  if (missing(ylab)) 
    ylab <- "Confidence interval for mean"
  plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, J + 0.3), xlab = xlab, 
       ylab = ylab, main = main)
  grid()
  points(ci.mat[,2], pch = pch, bg = "black")
  index <- 1:J
  segments(index, ci.mat[, 1], index, ci.mat[, 3], lwd = line.width[1])
  segments(index - 0.02, ci.mat[, 1], index + 0.02, ci.mat[, 1], lwd = line.width[2])
  segments(index - 0.02, ci.mat[, 3], index + 0.02, ci.mat[, 3], lwd = line.width[2])
  axis(1, 1:J, labels=level.labels)
}


plot.bsci.tiered <- function(data.frame, group.var=1, dv.var=2,  var.equal=FALSE, conf.level = 0.95, xlab = NULL, ylab = NULL, level.labels = NULL, main = NULL, pch = 19, pch.cex = 1.3, text.cex = 1.2, ylim = c(min.y, max.y), line.width= c(1.5, 1.5), tier.width=0, grid=TRUE) {
  data <- subset(data.frame, select=c(group.var, dv.var))
  fact <- factor(data[[1]])
  dv <- data[[2]]
  J <- nlevels(fact)
  ci.outer <- bsci(data.frame=data.frame , group.var=group.var, dv.var=dv.var,  difference=FALSE, var.equal=var.equal, conf.level =conf.level)
  ci.inner <- bsci(data.frame=data.frame , group.var=group.var, dv.var=dv.var,  difference=TRUE, var.equal=var.equal, conf.level =conf.level)
  moe.y <- max(ci.outer) - min(ci.outer)
  min.y <- min(ci.outer) - moe.y/3
  max.y <- max(ci.outer) + moe.y/3
  if (missing(xlab)) 
    xlab <- "Groups"
  if (missing(ylab)) 
    ylab <- "Confidence interval for mean"
  plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, J + 0.3), xlab = xlab, 
       ylab = ylab, main = main, cex.lab = text.cex)
  if (grid == TRUE) grid()
  points(ci.outer[,2], pch = pch, bg = "black", cex = pch.cex)
  index <- 1:J
  segments(index, ci.outer[, 1], index, ci.outer[, 3], lwd = line.width[1])
  axis(1, index, labels = level.labels)
  if(tier.width==0) {
    segments(index - 0.025, ci.inner[, 1], index + 0.025, ci.inner[, 1], lwd = line.width[2])
    segments(index - 0.025, ci.inner[, 3], index + 0.025, ci.inner[, 3], lwd = line.width[2])
  }
  else segments(index, ci.inner[, 1], index, ci.inner[, 3], lwd = line.width[1]*(1 + abs(tier.width)))
}


# meta-analysis of simple mean differences - functions from Baguley (submitted)
# These functions implement various forms of meta-analysis on simple (unstandardized) mean differences

md.ma.fixed <- function(mean.diffs, t.obs, df.err, conf.level = .95, NHST = FALSE){
  # fixed effects MA using Bond et al. (2003) method for simple mean differences
  std.errs <- mean.diffs/t.obs
  f.wts <- 1/(std.errs^2)
  mean.w <- sum((f.wts*mean.diffs))/sum(f.wts)
  J <- length(mean.diffs)
  l.wts <- vector("numeric",J)
  for(j in 1:J) l.wts[j] <- ((J-1)*f.wts[j] *(sum(f.wts)-f.wts[j]))/((J-1)*df.err[j]-4*(J-2))
  sigma.sq <- 1/sum(f.wts)*(1+4/sum(f.wts)^2*sum(l.wts))
  df <- sum(f.wts)^2/sum(f.wts^2/df.err)
  t.crit <- qt((1-conf.level)/2, df, lower.tail = FALSE)
  ci.mat <- matrix(,1,2, dimnames=list(paste(conf.level*100,'% CI'), c('lower', 'upper')))
  ci.mat[,1] <- mean.w - t.crit*sigma.sq^.5
  ci.mat[,2] <- mean.w + t.crit*sigma.sq^.5
  output <- list('number of studies in MA' = J, 'weighted mean difference' = mean.w, 'sampling variance' =sigma.sq, 'SE' = sigma.sq^.5, 'Interval' = ci.mat, 'df' = df)
  if(NHST==TRUE) output <- c(output, 't' = abs(mean.w/sigma.sq^.5), 'p' = pt(abs(mean.w/sigma.sq^.5), df, lower.tail = FALSE)*2)
  output
}


md.ma.random <- function(mean.diffs, t.obs, df.err, conf.level = .95, NHST = FALSE){
  # random effects MA using Hartung & Knapp (2001)/Bond et al. (2003) method for simple mean differences
  std.errs <- mean.diffs/t.obs
  f.wts <- 1/(std.errs^2)
  J <- length(mean.diffs)
  mean.uw = sum(mean.diffs)/J
  d.hat <- sum((f.wts*mean.diffs))/sum(f.wts)
  l.wts <- vector("numeric",J)
  for(j in 1:J)
    l.wts[j] <- ((J-1)*f.wts[j] *(sum(f.wts)-f.wts[j]))/((J-1)*df.err[j]-4*(J-2))
  d.hat.sigma.sq <- 1/sum(f.wts)*(1+4/sum(f.wts)^2*sum(l.wts))
  tau.sq <- var(mean.diffs)-(sum(1/f.wts))/J
  if(tau.sq < 0) tau.sq <- 0
  r.wts <- 1/(std.errs^2+tau.sq)
  d.tilda <- sum(r.wts*mean.diffs)/sum(r.wts)
  d.tilda.sigma.sq <- sum(r.wts*(mean.diffs-d.tilda)^2)/((J-1)*sum(r.wts))
  moe <- qt((1-conf.level)/2, J-1, lower.tail=FALSE) * d.tilda.sigma.sq^.5
  ci.mat <- matrix(,1,3, dimnames=list(paste(conf.level*100,'% CI'), c('lower', 'mean', 'upper')))
  ci.mat[,1:3] <- c(d.tilda-moe, d.tilda, d.tilda+moe)
  t.rand <- d.tilda/d.tilda.sigma.sq^.5
  p.rand <- pt(abs(t.rand), J-1, lower.tail=FALSE)*2
  output <- list('unweighted mean' = mean.uw, 'between study sampling variance' = tau.sq, 'total variance' = var(mean.diffs), 'R-square between' = tau.sq/var(mean.diffs), 'SE' = d.tilda.sigma.sq^.5, 'df' = J-1, 'interval estimate' = ci.mat)
  if(NHST==TRUE) output <- c(output, 't' = abs(t.rand), 'p' = p.rand)
  if(tau.sq <= 0) output <- c('Insufficient between study variation to estimate random effects', output)
  output
}


md.ma.het <- function(mean.diffs, t.obs, df.err) {
  # Bond et al. (2003) Fw test of heterogeneity of effect size in MA of simple mean differences
  std.errs <- mean.diffs/t.obs
  mean.uw = sum(mean.diffs)/length(mean.diffs)
  f.wts <- 1/(std.errs^2)
  mean.w <- sum(f.wts*mean.diffs)/sum(f.wts)
  J <- length(mean.diffs)
  u <- sum(1/df.err * (1-f.wts/sum(f.wts))^2)
  Fw <- (J+1)*sum(sum(f.wts*(mean.diffs-mean.w)^2))/(J^2-1+2*(J-2)*u)
  df1 <- J-1
  df2 <- (J^2-1)/(3*u)
  p.het <- pf(Fw, df1,df2, lower.tail = FALSE)
  output <- list('Heterogeneity test (Fw)' = Fw, 'df numerator' = df1, 'df denominator' = df2, 'p' = p.het)
  output
}

ma.t.contrast <- function(c.weights,  mean.diffs, t.obs, df.err, conf.level = .95, NHST = TRUE){
  # contrasts for Bond et al. (2003) MA of simple mean differences
  std.errs <- mean.diffs/t.obs
  f.wts <- 1/(std.errs^2)
  mean.w <- sum((f.wts*mean.diffs))/sum(f.wts)
  wj <- c.weights * 2/(sum(abs(c.weights)))
  t.contrast <- sum(wj*mean.diffs)/(sum(wj^2/f.wts)^.5)
  df.contrast <- sum(wj^2/f.wts)^2/sum(wj^4/(f.wts^2*df.err))
  moe <- qt((1-conf.level)/2, df.contrast, lower.tail=FALSE) * sum(wj^2/f.wts)^.5
  ci.mat <- matrix(,1,3, dimnames=list(paste(conf.level*100,'% CI'), c('lower', 'mean', 'upper')))
  ci.mat[,1:3] <- c(sum(wj*mean.diffs)-moe, sum(wj*mean.diffs), sum(wj*mean.diffs)+moe)
  p.contrast <- pt(abs(t.contrast), df.contrast, lower.tail=FALSE)*2
  r2.alerting <- cor(c.weights, mean.diffs)^2
  output <- list('contrast weights' = c.weights, 'interval estimate' = ci.mat, 'df' = df.contrast, 'r-square alerting' = r2.alerting)
  if(NHST==TRUE) output <- c(output, 't' = abs(t.contrast), 'p' = p.contrast)
  if(sum(c.weights)==0) output else 'contrast weights do not sum to zero'
}


paired.t.from.r <- function(mean1, mean2, sd1, sd2, n, r){
  # calculate paired t from summary data including correlation (r) between paired samples
  m.diffs <- mean1-mean2
  var.diffs <- sd1^2+sd2^2-2*r*sd1*sd2
  output <- m.diffs/((var.diffs)^.5/n^.5)
  output
}


md.vc.ma <- function(diff.ind, ind.sd1, ind.sd2, n1, n2 = n1, diff.paired, n.paired, t.paired, type = 'independent', conf.level = .95, NHST = FALSE){
  # Bonett (2009) varying coefficient MA of simple mean differences
  if(type == 'independent') {diff.paired = NULL ; n.paired = NULL ; t.paired = NULL}
  if(type == 'paired') {diff.ind = NULL ; ind.sd1 = NULL; ind.sd2 = NULL; n1 = NULL ; n2 = NULL}
  m <- length(diff.ind)+length(diff.paired)
  m.out <- m
  if(type == 'both') { m.out <- matrix(,1,2, dimnames = list('studies', c('independent', 'paired'))) ; m.out[1,] <- c(length(diff.ind), length(diff.paired))}
  phi.hat <- mean(c(diff.ind,diff.paired))
  var.phi <- c(ind.sd1^2/n1+ind.sd2^2/n2,(diff.paired/t.paired)^2)
  se.phi <- (m^-2*(sum(var.phi)))^.5
  numr.i <- sum(c(ind.sd1^2/n1,ind.sd2^2/n2))
  den.i <- sum(c(ind.sd1,ind.sd2)^4/((c(n1,n2))^3-(c(n1,n2))^2))
  numr.p <- sum((diff.paired/t.paired)^2)
  den.p <- sum(((((diff.paired/t.paired)^2)*n.paired)^2)/(n.paired^3-n.paired^2))
  df <- (numr.i+numr.p)^2/(den.i+den.p)
  moe <- se.phi*qt((1-conf.level)/2,df)*-1
  ci.mat <- matrix(,1,2, dimnames=list(paste(conf.level*100,'% CI'), c('lower', 'upper')))
  ci.mat[,1] <- phi.hat-moe
  ci.mat[,2] <- phi.hat+moe
  output <- list('number of studies in MA' = m.out, 'unweighted mean of differences' = phi.hat, 'interval' = ci.mat, 'df' = df)
  if(NHST==TRUE) output <- c(output, 't' = abs(phi.hat/se.phi), 'SE' = se.phi, 'p' = pt(abs(phi.hat/se.phi), df, lower.tail = FALSE)*2)
  output
}


md.vc.glm <- function(predictors, diff.ind, ind.sd1, ind.sd2, n1, n2 = n1, diff.paired, n.paired, t.paired, type = 'independent', conf.level = .95, NHST = FALSE){
  # GLM of heterogeneity of effects for Bonett (2009) varying coefficient MA of simple mean differences
  if(type == 'independent') {diff.paired = NULL ; n.paired = NULL ; t.paired = NULL}
  if(type == 'paired') {diff.ind = NULL ; ind.sd1 = NULL; ind.sd2 = NULL; n1 = NULL ; n2 = NULL}
  differences <- c(diff.ind, diff.paired)
  m.i <- length(diff.ind)
  m.p <- length(diff.paired)
  m <- length(differences)
  m.out <- m
  if(type == 'both') {m.out <- matrix(,1,2, dimnames = list('studies', c('independent', 'paired'))) ; m.out[1,] <- c(m.i, m.p)}
  if(is.vector(predictors)) {predictors <- as.data.frame(predictors) ; colnames(predictors) <- 'slope'}
  p.matrix <- as.matrix(predictors)
  p.names <- dimnames(p.matrix)
  full.names <- list(NULL, c('intercept', p.names[[2]]))
  q <- ncol(p.matrix)
  X.mat <- matrix(,m,q+1, dimnames = full.names)
  X.mat[,1] <- rep(1,m)
  for(i in 1:q) X.mat[,1+i] <- p.matrix[,i]
  beta.mat <- solve(t(X.mat) %*% X.mat) %*% (t(X.mat) %*% differences)
  colnames(beta.mat) <- 'coeffficient'
  var.phi <- c(ind.sd1^2/n1+ind.sd2^2/n2,(diff.paired/t.paired)^2)
  V.mat <- diag(var.phi)
  cov.mat <- solve(t(X.mat) %*% X.mat) %*% t(X.mat) %*% V.mat %*% X.mat %*% solve (t(X.mat)%*%X.mat)
  C.mat <- solve(t(X.mat) %*% X.mat) %*% t(X.mat)
  df.mat <- matrix(,1,q, dimnames = list(' ', p.names[[2]]))
  se.contrast <- matrix(,1,q)
  numr.i <- vector('numeric', length = q)
  den.i <- vector('numeric', length = q)
  numr.p <- vector('numeric', length = q)
  den.p <- vector('numeric', length = q)
  for(i in 1:q) numr.i[i] <- sum(C.mat[i+1, 1:m.i]^2*c(ind.sd1^2/n1, ind.sd2^2/n2))
  for(i in 1:q) den.i[i] <- sum(C.mat[i+1, 1:m.i]^4*c(ind.sd1,ind.sd2)^4/((c(n1,n2))^3-(c(n1,n2))^2))
  if(type != 'independent') {
    for(i in 1:q) numr.p[i] <- sum(C.mat[i+1, m.i+1:m.p]^2*(diff.paired/t.paired)^2);
    for(i in 1:q) den.p[i] <- sum(((((diff.paired/t.paired)^2)*n.paired)^2*C.mat[i+1, m.i+1:m.p]^4)/(n.paired^3-n.paired^2))
  }
  for(i in 1:q) df.mat[i] <- (numr.i[i]+numr.p[i])^2/(den.i[i]+den.p[i])
  for(i in 1:q) se.contrast[i] <- (sum(C.mat[i+1,]^2*var.phi))^.5
  t.crit <- qt((1-conf.level)/2,df.mat[1,])*-1
  ci.mat <- matrix(,2,q, dimnames=list(c(paste(conf.level*100,'% CI lower'), paste(conf.level*100,'% CI upper')), p.names[[2]]))
  for(i in 1:q) {ci.mat[2,i] <- beta.mat[i+1,] + t.crit[i]*diag(cov.mat)[i+1]^.5;  ci.mat[1,i] <- beta.mat[i+1,] - t.crit[i]*diag(cov.mat)[i+1]^.5}
  t.mat <- matrix(,1,q+1)
  p.mat <- matrix(,1,q)
  t.mat <- abs(t(beta.mat)/diag(cov.mat)^.5)
  t.mat <- t(as.matrix(t.mat[,1:q+1]))
  dimnames(t.mat) <- list('t', p.names[[2]])
  p.mat <- pt(t.mat, df.mat, lower.tail = FALSE)*2
  dimnames(p.mat) <- list('p', p.names[[2]])
  output <- list('number of studies in MA' = m.out, 'Parameter estimates' = beta.mat, 'Asymptotic covariance matrix' = cov.mat, 'intervals' = ci.mat, 'df' = df.mat, 'standard errors' = diag(cov.mat)^.5)
  if(q==1) output <- list('number of studies in MA' = m.out, 'Parameter estimates' = beta.mat, 'Asymptotic covariance matrix' = cov.mat, 'intervals' = ci.mat, 'df' = df.mat, 'standard errors' = diag(cov.mat)^.5, 'r-square alerting' = cor(p.matrix[,1], differences)^2)
  if(NHST==TRUE) output <- list('number of studies in MA' = m.out, 'Parameter estimates' = beta.mat, 'Asymptotic covariance matrix' = cov.mat, 'intervals' = ci.mat, 'df' = df.mat, 'standard errors' = diag(cov.mat)^.5, 'test statistic' = t.mat, 'statistical significance' = p.mat)
  output
}


decom.plot <- function(mean.diffs, sd.diffs, xlab, ylab, main, leg.x, leg.y, leg.cex = .85, mar = c(7,4,3,2)){
  # decomposition plot after Bond et al. (2003)
  par(mar=mar)
  g.slope <- mean(mean.diffs/sd.diffs)
  range.m <- range(mean.diffs)[2]-range(mean.diffs)[1]
  if (missing(main)) main <- 'Decomposition plot'
  if (missing(ylab)) ylab <- 'Simple (raw) mean differences (G)'
  if (missing(xlab)) xlab <- 'Standard deviation of the differences'
  if (missing(leg.x)) leg.x <- mean(sd.diffs)*.8
  if (missing(leg.y)) leg.y <- min(mean.diffs)-range.m/4
  plot(sd.diffs, mean.diffs, xlab = xlab, ylab=ylab, main=main, asp = 1, cex.main = 1)
  abline(mean(mean.diffs), 0, lty=2)
  abline(0, g.slope, lty=3)
  legend(leg.x, leg.y, legend = c('Unweighted mean G', 'Unweighted mean g'), lty=c(2,3), cex= leg.cex, bty='n', xpd = TRUE)
}

md.vc.resid <- function(predictors, diff.ind =NULL, ind.sd1=NULL, ind.sd2=NULL, n1=NULL, n2 = n1, diff.paired =NULL, n.paired=NULL, t.paired=NULL, type = 'independent', plot = TRUE, reps =2000, bonferroni = FALSE, conf.level = .95, NHST = FALSE, pch = 1){
  # extracts residuals, standardized residuals, leverage and normal probability plot for varying coefficient MA of simple mean differences
  diffs <- c(diff.ind,diff.paired)
  m <- length(diffs)
  if(missing(predictors)) {beta.mat <- matrix(mean(diffs)) ; X.mat <-  matrix(1, nrow = m, ncol=1) } else {
    beta.mat <- rfe.glm(predictors = predictors, diff.ind = diff.ind, ind.sd1=ind.sd1, ind.sd2=ind.sd2, n1=n1, n2 = n2, diff.paired = diff.paired, n.paired = n.paired, t.paired = t.paired, type = type)[[2]] ;
    if(is.vector(predictors)) {predictors <- as.data.frame(predictors) ; colnames(predictors) <- 'slope'}
    p.matrix <- as.matrix(predictors)
    q <- ncol(p.matrix)
    X.mat <- matrix(,m,q+1)
    X.mat[,1] <- rep(1,m)
    for(i in 1:q) X.mat[,1+i] <- p.matrix[,i]
  }
  h.mat <- X.mat %*% solve(t(X.mat) %*% X.mat) %*% t(X.mat)
  var.phi <- c(ind.sd1^2/n1+ind.sd2^2/n2,(diff.paired/t.paired)^2)
  v.mat <- diag(var.phi)
  # portions of the following R code are adapted from the Viechtbauer (2010) metafor package qqnorm.rma.uni() function
  ve.mat <- (diag(m)-h.mat) %*% tcrossprod(v.mat, (diag(m)-h.mat))
  se.resid <- sqrt(diag(ve.mat))
  y.hat <- X.mat %*% beta.mat
  resids <- as.vector(diffs-y.hat)
  output <- list('residuals' = resids, 'predicted' = as.vector(y.hat), 'observed' = diffs, 'z' = resids/se.resid, 'leverage' = diag(h.mat))
  if(plot) {qqnorm(resids/se.resid, pch =pch) ; abline(0,1);
    dat <- matrix(rnorm(m * reps), nrow = m, ncol = reps)
    alpha <- 1-conf.level;
    ei <- (diag(m)-h.mat) %*% dat;
    ei <- apply(ei, 2, sort);
    if (bonferroni) {
      lb <- apply(ei, 1, quantile, (alpha/2)/m)
      ub <- apply(ei, 1, quantile, 1 - (alpha/2)/m)
    }
    else {
      lb <- apply(ei, 1, quantile, (alpha/2))
      ub <- apply(ei, 1, quantile, 1 - (alpha/2))
    }
    temp <- qqnorm(lb, plot.it = FALSE)
    temp <- supsmu(temp$x, temp$y)
    lines(temp$x, temp$y, lty = 3)
    temp <- qqnorm(ub, plot.it = FALSE)
    temp <- supsmu(temp$x, temp$y)
    lines(temp$x, temp$y, lty = 3)
  }
  output
}


# mac OS paste data from clipboard (tab limited columns) - see URL <to be added>

paste.data <- function(header=FALSE) {read.table(pipe("pbpaste"), header=header)}

# additional functions for prior exposure workshops 2015

# smd.post() calculates a posterior distribution for a standardized mean difference
# at present the function just works for smd from 2 independent means
# input is the observed smd (assumed to be Hedges' g) and n per group


smd.post <- function(smd.obs, smd.prior, sigma.prior, n1, n2=n1, probability = 0.95, plot=FALSE){
  sigma.obs <- ((n1 + n2)/(n1 * n2) + smd.obs^2/(2 * (n1 + n2 - 2)))^0.5
  precision.prior <- 1/sigma.prior^2
  precision.lik <- 1/sigma.obs^2
  precision.post <- precision.prior + precision.lik
  smd.post <- (precision.prior/precision.post) * smd.prior + (precision.lik/precision.post) * smd.obs
  sigma.post <- (1/precision.post)^0.5
  lower <- smd.post + qnorm((1 - probability)/2) * sigma.post
  upper <- smd.post - qnorm((1 - probability)/2) * sigma.post
  post.prob.int <- c(lower, upper)
  min.m <- min(smd.obs, smd.prior, smd.post)
  max.m <- max(smd.obs, smd.prior, smd.post)
  max.s <- max(sigma.obs, sigma.prior, sigma.post)
  y.max <- max(dnorm(smd.prior, smd.prior, sigma.prior), dnorm(smd.obs, smd.obs, sigma.obs), dnorm(smd.post, smd.post, sigma.post))
  if (plot == TRUE) {
    curve(dnorm(x, smd.post, sigma.post), lty = 3, lwd=3, xlim = c(min.m - 3.75 * max.s, max.m + 3.75 * max.s), yaxt = "n", ylab = NA, xlab = expression(delta), col='purple')
    curve(dnorm(x, smd.prior, sigma.prior), lwd = 2, col='dark blue', add = TRUE)
    curve(dnorm(x, smd.obs, sigma.obs), lwd=2, col='orange', add = TRUE)
    legend('topright', y.max, legend = c("Likelihood", "Prior", "Posterior"), lty = c(1, 1, 3), lwd = c(2, 2, 3), col=c('orange', 'dark blue', 'purple'))
  }
  output <- list(`Posterior smd` = smd.post, `Posterior SD` = sigma.post, `Probability interval` = post.prob.int, Probability = probability, lower = lower, upper = upper)
  output
}


