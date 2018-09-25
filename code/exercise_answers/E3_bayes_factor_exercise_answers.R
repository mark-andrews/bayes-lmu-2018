

# Bayes factor exercise
# note needs functions from Serious Stats book and BayesFactor package
source('http://www2.ntupsychology.net/seriousstats/SeriousStatsAllfunctions.txt')
library(BayesFactor)

##############################
# Bayes factor exercise Q1a
##############################

unit.prior.Bf.2s(t=2.02, n1=91.5, n2=91.5, scale.factor=1)


##############################
# Bayes factor exercise Q1b 
##############################

unit.prior.Bf.2s(t=2.02, n1=91.5, n2=91.5, scale.factor=0.5)

unit.prior.Bf.2s(t=2.02, n1=91.5, n2=91.5, scale.factor=2)

##############################
# Bayes factor exercise Q1c
##############################

unit.prior.Bf.2s(t=2.02, n1=91.5, n2=91.5, scale.factor=0.05)

unit.prior.Bf.2s(t=2.02, n1=91.5, n2=91.5, scale.factor=0.01)

unit.prior.Bf.2s(t=2.02, n1=91.5, n2=91.5, scale.factor=0.001)

unit.prior.Bf.2s(t=2.02, n1=91.5, n2=91.5, scale.factor=100)

# this illustrates the automatic Occam's razor property of Bayes factors
# as you become more certain of your prior (its precision is greater/variance smaller)
# the posterior depends on whether the prior matches the data
# so a strong certain prior may decrease the BF if the data contradict it


##############################
# Bayes factor exercise Q2a
##############################

exp(ttest.tstat(t=2.02, n1=91.5, n2=91.5, rscale='medium')[['bf']])

# the log is used internally within the BayesFactor functions in order to get accurate calculations
# with very large or small numbers


##############################
# Bayes factor exercise Q2b
##############################

exp(ttest.tstat(t=2.02, n1=91.5, n2=91.5, rscale='medium', nullInterval=c(0,Inf))[['bf']])



##############################
# Bayes factor exercise Q2c
##############################

exp(ttest.tstat(t=2.02, n1=91.5, n2=91.5, rscale='medium', nullInterval=c(-Inf,0))[['bf']])

##############################
# Bayes factor exercise Q2d
##############################

exp(ttest.tstat(t=2.02, n1=91.5, n2=91.5, rscale='medium', nullInterval=c(-Inf,Inf))[['bf']])



##############################
# Bayes factor exercise Q2e
##############################

bf1 <- exp(ttest.tstat(t=2.02, n1=91.5, n2=91.5, rscale='medium', nullInterval=c(-Inf,0))[['bf']])
bf2 <- exp(ttest.tstat(t=2.02, n1=91.5, n2=91.5, rscale='medium', nullInterval=c(0,Inf))[['bf']])

bf2/bf1

# this is a test of whether the effect is positive versus negative (and assumes implicitly that the effect is non-zero)

##############################
# Bayes factor exercise OPTIONAL
##############################

# some examples here using raw data
# these commands load up a simple data set (from a experiment with four experimental conditions and some covariates)

diag.data <- read.csv("http://www2.ntupsychology.net/seriousstats/diagram.csv")

# this extracts two conditions into vectors to use in a t test
cont <- diag.data$descript[1:10]
seg <- diag.data$descript[31:40]

# frequentist t test

t.test(seg, cont)

# BayesFactor t test with default JZS priors
ttestBF(seg, cont)

# directional (one-tailed) version of Bayesian t test
ttestBF(seg, cont, nullInterval=c(0, Inf))


# Bayesian ANOVA with one factor
bf.anov <- anovaBF(descript ~ group, diag.data)
bf.anov

# Bayes Factor simple linear regression 
regressionBF(descript ~ time, data = diag.data, rscaleCont = 'medium')

# ANCOVA
bf.ancov <- lmBF(descript ~ time + group, data = diag.data)
bf.ancov


# brms example

library(brms)

anov.brm1 <- brm(descript ~ 1, diag.data)
anov.brm1
plot(anov.brm1)

anov.brm2 <- brm(descript ~ group, diag.data, sample_prior=TRUE)
anov.brm2
plot(anov.brm2)

waic(anov.brm1, anov.brm2)
loo(anov.brm1, anov.brm2)


hypothesis(anov.brm2, "groupPicture < groupSegmentedDiagram")
hypothesis(anov.brm2, "groupPicture > groupText")
hypothesis(anov.brm2, "groupSegmentedDiagram + Intercept > Intercept")
hypothesis(anov.brm2, "groupSegmentedDiagram + Intercept > Intercept + groupText + groupPicture")
hypothesis(anov.brm2, "groupSegmentedDiagram  > groupText + groupPicture")



