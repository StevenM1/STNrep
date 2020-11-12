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


# some descriptives
dat$acc <- as.numeric(dat$R) == as.numeric(dat$S)
RTbyCoh <- aggregate(RT~s*coherence, dat, mean)
accByCoh <- aggregate(acc~s*coherence, dat, mean)

par(mfrow=c(1,2))
boxplot(RT~coherence, RTbyCoh)
points(factor(RTbyCoh$coherence), RTbyCoh$RT, col=1)
boxplot(acc~coherence, accByCoh)
points(factor(accByCoh$coherence), accByCoh$acc, col=1)


# by ds
dat$ds <- 1
dat$ds[as.numeric(dat$s) > 19] <- 2
dat$ds <- as.factor(dat$ds)
RTbyCoh <- aggregate(RT~s*coherence*ds, dat, mean)
accByCoh <- aggregate(acc~s*coherence*ds, dat, mean)
par(mfrow=c(1,2))
boxplot(RT~coherence*ds, RTbyCoh)
boxplot(acc~coherence*ds, accByCoh)

# differences?
library(reshape2)
RTbyCohWide <- reshape(RTbyCoh, direction='wide', v.names='RT', timevar='coherence', idvar=c('s', 'ds'))
RTbyCohWide$diff <- RTbyCohWide$RT.low-RTbyCohWide$RT.high
RTbyCohLong <- melt(RTbyCohWide, value.name='RT')
RTbyCohLong <- RTbyCohLong[RTbyCohLong$variable == 'diff',]
boxplot(RT~ds, RTbyCohLong)
points(factor(RTbyCohLong$ds), RTbyCohLong$RT, col=1)

accByCohWide <- reshape(accByCoh, direction='wide', v.names='acc', timevar='coherence', idvar=c('s', 'ds'))
accByCohWide$diff <- accByCohWide$acc.high-accByCohWide$acc.low
accByCohLong <- melt(accByCohWide, value.name='acc')
accByCohLong <- accByCohLong[accByCohLong$variable == 'diff',]
boxplot(acc~ds, accByCohLong)
points(factor(accByCohLong$ds), accByCohLong$acc, col=1)

