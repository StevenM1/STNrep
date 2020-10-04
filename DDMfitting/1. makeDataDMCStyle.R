# load data
dat <- read.csv('../derivatives/behavior.csv')
dat[dat$ds=='ds-02','subject'] = dat[dat$ds=='ds-02','subject'] + 19
dat$s <- factor(dat$subject)
dat$S <- factor(ifelse(dat$stimulus == 'left', 'stimleft', 'stimright'))
dat$coherence[dat$coherence == 0.16] = 'high'
dat$coherence[dat$coherence == 0.35] = 'high'
dat$coherence[dat$coherence == 0.08] = 'low'
dat$coherence[dat$coherence == 0.15] = 'low'
dat$coherence <- factor(dat$coherence, levels=c('low', 'high'))
dat$cue <- paste0('cue', dat$cue)
dat$cue <- factor(dat$cue, levels=c('cueleft', 'cueright', 'cueneutral'))
dat$RT <- dat$rt/1000
dat <- dat[!dat$response == -1,]
dat$R <- factor(ifelse(dat$response == 1, 'r1', 'r2'), levels=c('r1', 'r2'))
dat <- dat[,c('s', 'S', 'R', 'RT', 'coherence', 'cue')]
# remove extremely fast RTs
dat <- dat[dat$RT>.15,]
save(dat, file='behavior.RData')