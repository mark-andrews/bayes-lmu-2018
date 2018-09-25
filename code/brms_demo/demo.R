library(brms)
library(dplyr)
library(readr)

Df <- read_csv("../../data/diagram.csv") %>% 
  select(-X1)

M_null <- brm(descript ~ 1, Df, save_all_pars = T)
summary(M_null)
plot(M_null)

M <- brm(descript ~ group, Df, save_all_pars = T)
summary(M)

plot(M)

# Out of sample predictive performance 
waic(M_null, M)
loo(M_null, M)

# Bayes factor of M_null v M
bayes_factor(M, M_null)

# Directional hypothesis testing
hypothesis(M, "groupPicture < groupSegmentedDiagram")
hypothesis(M, "groupPicture > groupText")
hypothesis(M, "groupSegmentedDiagram + Intercept > Intercept")
hypothesis(M, "groupSegmentedDiagram + Intercept > Intercept + groupText + groupPicture")
hypothesis(M, "groupSegmentedDiagram  > groupText + groupPicture")