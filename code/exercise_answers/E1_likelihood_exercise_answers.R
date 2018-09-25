

# simple likelihood examples


# these use R distribution functions
###############################################
# OPTIONAL intro on distribution functions in R
###############################################

# get help on distribution functions
# what distributions are available in base R (there are packages that expand this)

?distributions

# binomial distribution functions
?dbinom

# using distribution functions
# here we get the probability of getting exactly 19 out 25 2AFC responses correct assuming P =.5 (guessing)
# d here is for density - but note that for discrete distributions this is a probability mass not probability density

dbinom(19, 25, .5)

# cumulative probability - lower or upper tail

pbinom(18, 25, .5, lower.tail=FALSE)
pbinom(18, 25, .5, lower.tail=TRUE)

# this is equivalent to summing probability density/mass up to this point
# so the following quantities are identical

1 - dbinom(c(19:25), 25, .5)
pbinom(18, 25, .5, lower.tail=FALSE)

# we can also easily get quantiles
# for a discrete function this returns the highest value that includes at least this proportion of the lower tail (or lowest value for upper tail)
qbinom(.975, 25, .5)

# each function also has a random generation function (useful for simulations etc.) 

rbinom(10, 1, .5)
rnorm(n = 10, mean = 100, sd = 15)

# all these functions can be easily plotted e.g., using curve() or plot(..., type='h')

curve(dnorm(x, 100, 15), xlim = c(45, 155), ylab = "Probability density") # 'x' mandatory arg for curve()

down vo
plot(0:10, dbinom(0:10, 10, .5), type="h", xlab = 'Number of heads')

# fancier example
plot(0:20, dpois(0:20, 6), type="h", main="Poisson probability mass function with approximating normal density")
abline(h=0, col="green2")
curve(dnorm(x, 6, sqrt(6)), lwd=2, col="red", add=TRUE)

##############################
# Likelihood exercise EXAMPLE
##############################
# Likelihood for 19 successes on 25 trials of a 2AFC task for observed probability

.76^19 * .24^6

# Likelihood for 19 successes on 25 trials of a 2AFC task if guessing
.5^19 * .5^6

# Likelihood ratio of guessing versus observed performance

(.76^19 * .24^6)/(.5^19 * .5^6)

# example repeated using distribution functions
# note that the likelihoods have changed but the ratio has not
# a likelihood is probability of parameter given data or any quantity proportional to it
# as long as we use the same constant to scale the probabilities the ratios will be correct

dbinom(19, 25, .76) # note that this is the maximum likelihood estimate
dbinom(19, 25, .50)
dbinom(19, 25, .76) / dbinom(19, 25, .50)

##############################
# Likelihood exercise Q1
##############################
# likelihood of typical performance versus observed

dbinom(19, 25, .80)
dbinom(19, 25, .80) / dbinom(19, 25, .76)

##############################
# Likelihood exercise Q2
##############################
# LR of typical performance versus guessing
# this is a much better way to test a likelihood ratio as it doesn't rely on the MLE (which maximizes one of the inputs)
dbinom(19, 25, .80) / dbinom(19, 25, .50)

##############################
# Likelihood exercise Q3
##############################

# the only term that varies is the third term (the probability)
# this is the unknown parameter value
# likelihood treats the data as known/observed and the parameter as unknown
# hence we are comparing hypotheses about parameters given the observed data
# this is at the heart of the evidential / likelihood approach to inference


##############################
# Likelihood exercise Q4
##############################

# you can use any data here - perhaps contrast with an analysis that uses a t test
# this example uses data from the flag priming study
LR.t(mu1 = 0, mu2 = 0.06, mu.obs = 0.142, SE = 0.0703, df = 183)

##############################
# Likelihood exercise Q5
##############################

# here is an opportunity to try out a function to obtain and plot likelihood intervals

t.lik.int(mean = 10, SE = 2, df = 28, independent = TRUE, plot=TRUE)


##############################
# Optional stopping
##############################

# this produces variations of my optional stopping plot for one sample t when H0 = TRUE
# the green line plots the likelihood ratio which increases gradually to favour the null

# this is an initial sample of n = 10 and adding 1 data point before calculating a p value and likelihood ratio 
t.updated(initial.n=10, increment=1, iterations=1000, plot=TRUE)

# this is an initial sample of n = 30 and adding 1 data point before calculating a p value and likelihood ratio 
t.updated(initial.n=100, increment=5, iterations=1000, plot=TRUE)
# note that increasing starting n helps delay the impact of optional stopping (but p is still biased)










