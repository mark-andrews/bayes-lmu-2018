alpha <- 1
beta <- 1
m <- 139
n <- 250

# Log marginal likelihood of Model = Beta(alpha, beta) when n = 250, m = 139.
log.marginal.likelihood.model <- lbeta(alpha + m, n - m + beta)
log.marginal.likelihood.alt.model <- lgamma(alpha + m) + lgamma(beta + n - m) - lgamma(alpha + beta + n)
log.marginal.likelihood.alt.2.model <- lfactorial(m) + lfactorial(n-m) - lfactorial(n + 1)

# Log marginal likelihood of Null Model (theta = 1/2) when n = 250, m = 139.
log.marginal.likelihood.null <- -(n * log(2))

# Log Bayes factor of Model = Beta(alpha, beta) to Model = (theta = 1/2).
log.bayes.factor <- log.marginal.likelihood.model - log.marginal.likelihood.null
bayes.factor <- exp(log.bayes.factor )

# Relative probability of null model
exp(log.marginal.likelihood.null - log.marginal.likelihood.model)

classical.test <- binom.test(x=m, n=n, p = 0.5)