library(rstanarm)
library(ggplot2)

load('../../data/beautyeval.Rda')

# Visualize it
ggplot(beautydata,
       mapping=aes(x=beauty, y=eval, col=sex)) + 
  geom_point() +
  geom_smooth(method='lm') +
  xlab('Lecturer attractiveness') +
  ylab('Teaching evaluation score') +
  ggtitle('Do good looking lecturers get better teaching evaluation score?')


# Classical glm
M.lm <- lm(eval ~ beauty*sex, data=beautydata)

# Linear model with Stan 
M.stan.lm <- stan_lm(eval ~ beauty*sex, 
                  prior = R2(location=0.5),
                  chains = 3,
                  cores = 3,
                  seed = 1001,
                  data=beautydata)

# Linear model with Stan GLM
M.stan.glm <- stan_glm(eval ~ beauty*sex,
                   data = beautydata, 
                   family = gaussian(), 
                   prior = cauchy(), 
                   prior_intercept = cauchy(),
                   chains = 3, 
                   cores = 3, 
                   seed = 42)

# Linear model with Stan GLM
M.stan.glm.additive <- stan_glm(eval ~ beauty + sex,
                                data = beautydata, 
                                family = gaussian(), 
                                prior = cauchy(), 
                                prior_intercept = cauchy(),
                                chains = 3, 
                                cores = 3, 
                                seed = 101)

# Compare additive and full models 
compare_models(waic(M.stan.glm), waic(M.stan.glm.additive))

# Posterior predictive checks
pp_check(M.stan.glm.additive, plotfun = "hist", nreps = 8, binwidth=0.1)
pp_check(M.stan.glm.additive, plotfun = "stat", stat = "mean")
pp_check(M.stan.glm.additive, plotfun = "stat_2d", stat = c("mean", "sd"))
