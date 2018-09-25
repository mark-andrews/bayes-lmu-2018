library(lme4)
library(rjags)
library(MCMCglmm)

pitch.df <- read.csv('pitch.csv', header=T)

participant <- pitch.df$Participant
pitch <- pitch.df$pitch
base <- pitch.df$base
attract <- pitch.df$attract
face <- pitch.df$Face

N <- length(participant)
J <- max(participant)
K <- max(face)

n.iter <- 50000

### Model 1: pitch ~ 1 + (1|Participant) ###

M.1 <- jags.model("model.1.jags", 
                  data = list('pitch'=pitch, 'J'=J, 'N'=N, 'participant'=participant),
                  n.chain=3)

update(M.1, n.iter)
S.1 <- coda.samples(M.1, c('alpha.0', 'sigma', 'sigma.0'), n.iter)

M.1.lmer <- lmer(pitch ~ 1 + (1|Participant), pitch.df)
pitch.mcmc.1 <- MCMCglmm(pitch ~ 1, 
                         random= ~ Participant, 
                         nitt=n.iter, 
                         data=pitch.df)

summary(S.1)
summary(M.1.lmer)
summary(pitch.mcmc.1)

### Model 2: pitch ~ base + (1|Participant) ###

M.2 <- jags.model("model.2.jags", 
                  data = list('pitch'=pitch, 'J'=J, 'N'=N, 'participant'=participant, 'base'=base),
                  n.chain=3)
update(M.2, n.iter)
S.2 <- coda.samples(M.2, c('alpha.0', 'beta.base', 'sigma', 'sigma.0'), n.iter)

M.2.lmer <- lmer(pitch ~ base + (1|Participant), pitch.df)

pitch.mcmc.2 <- MCMCglmm(pitch ~ base, 
                         random= ~ Participant, 
                         nitt=n.iter, 
                         data=pitch.df)

summary(S.2)
summary(M.2.lmer)
summary(pitch.mcmc.2)

### Model 3: pitch ~ base + attract + (1|Participant) ###

M.3 <- jags.model("model.3.jags", 
                  data = list('pitch'=pitch, 'J'=J, 'N'=N, 'participant'=participant, 'base'=base, 'attract'=attract),
                  n.chain=3)

update(M.3, n.iter)
S.3 <- coda.samples(M.3, c('alpha.0', 'beta.base', 'beta.attract', 'sigma', 'sigma.0'), n.iter)

M.3.lmer <- lmer(pitch ~ base + attract + (1|Participant), pitch.df)
pitch.mcmc.3 <- MCMCglmm(pitch ~ base + attract, 
                         random= ~ Participant, 
                         nitt=n.iter, 
                         data=pitch.df)

summary(S.3)
summary(M.3.lmer)
summary(pitch.mcmc.3)

### Model 4: pitch ~ base + attract + (attract|Participant) ###

M.4 <- jags.model("model.4.jags", 
                  data = list('pitch'=pitch, 'J'=J, 'N'=N, 'participant'=participant, 'base'=base, 'attract'=attract),
                  n.chain=3)

update(M.4, n.iter)
S.4 <- coda.samples(M.4, 
                    c('alpha.0', 'beta.base', 'beta.attract.0', 'sigma', 'sigma.alpha.0', 'sigma.beta.attract.0'), 
                    n.iter)

M.4.lmer <- lmer(pitch ~ base + attract + (attract|Participant), pitch.df)
pitch.mcmc.4 <- MCMCglmm(pitch ~ base + attract, 
                         random= ~ idh(1+attract):Participant, 
                         nitt=n.iter, 
                         data=pitch.df)


summary(S.4)
summary(M.4.lmer)
summary(pitch.mcmc.4)

### Model 5: pitch ~ base + attract + (1|Participant) + (1|Face) ###
M.5 <- jags.model("model.5.jags", 
                  data = list('pitch'=pitch, 'J'=J, 'N'=N, 'K'=K, 'participant'=participant, 'base'=base, 'attract'=attract, 'face'=face),
                  n.chain=3)
update(M.5, n.iter)
S.5 <- coda.samples(M.5, 
                    c('alpha',
                      'beta.base', 
                      'beta.attract', 
                      'sigma',
                      'sigma.face.alpha.0', 
                      'sigma.participant.alpha.0'), 
                    n.iter)

M.5.lmer <- lmer(pitch ~ base + attract + (1|Participant) + (1|Face), data=pitch.df)
pitch.mcmc.5 <- MCMCglmm(pitch ~ base + attract, 
                         random= ~ Participant + Face, 
                         nitt=n.iter, 
                         data=pitch.df)

### Model 6: pitch ~ base + attract + (attract|Participant) + (1|Face) ###
M.6 <- jags.model("model.6.jags", 
                  data=list('pitch'=pitch, 'J'=J, 'N'=N, 'K'=K, 'participant'=participant, 'base'=base, 'attract'=attract, 'face'=face),
                  n.chain=3)
update(M.6, n.iter)
S.6 <- coda.samples(M.6, 
                    c('alpha',
                      'beta.base', 
                      'beta.attract.0', 
                      'sigma',
                      'sigma.face.alpha.0', 
                      'sigma.participant.alpha.0'), 
                    n.iter)

M.6.lmer <- lmer(pitch ~ base + attract + (attract|Participant) + (1|Face), data=pitch.df)
pitch.mcmc.6 <- MCMCglmm(pitch ~ base + attract, 
                         random= ~ idh(1+attract):Participant + Face, 
                         nitt=n.iter, 
                         data=pitch.df)


### Model 7: pitch ~ base + attract + (base + attract|Participant) + (1|Face) ###
M.7 <- jags.model("model.7.jags", 
                  data=list('pitch'=pitch, 'J'=J, 'N'=N, 'K'=K, 'participant'=participant, 'base'=base, 'attract'=attract, 'face'=face),
                  n.chain=3)
update(M.7, n.iter)
S.7 <- coda.samples(M.7, 
                    c('alpha',
                      'beta.base.0', 
                      'beta.attract.0', 
                      'sigma',
                      'sigma.face.alpha.0', 
                      'sigma.participant.alpha.0',
                      'sigma.beta.base.0',
                      'sigma.beta.attract.0'), 
                    n.iter)

M.7.lmer <- lmer(pitch ~ base + attract + (attract+base|Participant) + (1|Face), data=pitch.df)
M.7.lmer <- lmer(pitch ~ base + attract 
                 + (1|Participant) 
                 + (0+ base|Participant) 
                 + (0+ attract|Participant) 
                 + (0 + attract : base|Participant) 
                 + (1|Face), data=pitch.df)

pitch.mcmc.7 <- MCMCglmm(pitch ~ base + attract, 
                         random= ~ idh(1 + base+attract):Participant + Face, 
                         nitt=n.iter, 
                         data=pitch.df)
