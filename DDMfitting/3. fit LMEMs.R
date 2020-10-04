load('behavior.RData')
library(lme4)
library(lmerTest)
library(report)

head(dat)
dat$cue_congruency <- 'Neutral'
dat[dat$S=='stimleft'&dat$cue=='cueleft','cue_congruency'] <- 'Congruent'
dat[dat$S=='stimright'&dat$cue=='cueright','cue_congruency'] <- 'Congruent'
dat[dat$S=='stimleft'&dat$cue=='cueright','cue_congruency'] <- 'Incongruent'
dat[dat$S=='stimright'&dat$cue=='cueleft','cue_congruency'] <- 'Incongruent'
dat$cue_congruency <- factor(dat$cue_congruency, levels=c('Neutral', 'Congruent', 'Incongruent'))
dat$accuracy <- factor(as.numeric(dat$R)==as.numeric(dat$S), levels = c(TRUE, FALSE))

lmer1 <- lmer(RT~coherence*cue_congruency + (1|s), dat)
summary(lmer1)

lmer1 <- lmer(RT~coherence*cue_congruency*accuracy + (1|s), dat)
summary(lmer1)
r <- report(lmer1)
table_long(r)

dat$accuracy <- factor(as.numeric(dat$R)==as.numeric(dat$S))#, levels = c(TRUE, FALSE))
lmer2 <- glmer(accuracy ~ coherence*cue_congruency + (1|s), dat, family='binomial')
summary(lmer2)
r <- report(lmer2)
table_long(r)

dat$accuracy <- as.logical(dat$accuracy)
aggregate(accuracy~s, dat, mean)
