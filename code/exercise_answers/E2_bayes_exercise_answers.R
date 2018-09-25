

# Bayes exercise examples


##############################
# Bayes exercise Q1
##############################
# here I have set up the calculations - note that the likelihood SD is the SE (standard error) from the data
# this follows from the definition of a standard error as the SD of a sampling distribution for a paramater such as a mean
# so the SE is in a general sense the SD of a parameter

# setting up the mean and variance of the likelihood and the prior distributions
m.obs <- 9.6
var.obs <- 4.7^2
m.prior <- 2.0
var.prior <- 5.0^2

# setting up the calculations
# not all Bayesian analyses have such simple analytic results
# Also many complex models don't have (known) analytic solutions

var.post <- 1/(1/var.obs + 1/var.prior)
m.post <-  var.post/var.obs * m.obs + var.post/var.prior * m.prior
m.post ; var.post

# 95% posterior probability (credibility) interval 
m.post + c(-1,1) * qnorm(.975) * sqrt(var.post)



##############################
# Bayes exercise Q2
##############################

# With the added knowledge that N = 30 one can calculate the sample SD from the SE

s.obs <- sqrt(var.obs * 30)
s.obs

Bayes.norm.1s(9.6, 25.74, 2, 5, N = 30, ssc=FALSE, plot=TRUE)

# with small sample correction
Bayes.norm.1s(9.6, 25.74, 2, 5, N = 30, ssc=TRUE, plot=TRUE) 

##############################
# Bayes exercise OPTIONAL
##############################

# try out the two sample test - this just uses the mean differences
# Bayes.norm.2s(M.diff, SD.1, SD.2, diff.prior, sigma.prior, n1, n2, probability = 0.95, ssc = TRUE, plot = FALSE)
# for example

Bayes.norm.2s(5, 3, 4, diff.prior=10, sigma.prior=10, 32, 11, probability = 0.95, ssc = TRUE, plot = FALSE)

