rm(list=ls())
setwd('~/bias/DDMfitting')
source ("dmc/dmc.R")
load_model ("DDM","ddm.R") 
source('./model3.R')
source('./plotUtils.R')

## Load winning model
load('samplesModel3.RData')
load('behavior.RData')

# Simulate model
ppSim <- h.post.predict.dmc(samples, cores=30, save.simulation = TRUE)

# Get data, make accuracy column
data <- lapply(samples, function(x) x$data)
ppSim <- lapply(ppSim, getAccuracy)
data <- lapply(data, getAccuracy)

# Quantiles in model predictions
qps <- seq(.1, .9, .2)
pp <- lapply(ppSim, getQuantiles, qps=qps)
qpsCorrect <- do.call(rbind, lapply(pp, function(x) x$qpsCorrect))
qpsError <- do.call(rbind, lapply(pp, function(x) x$qpsError))
accuracies <- do.call(rbind, lapply(pp, function(x) x$accuracy))
qpsByRep <- do.call(rbind, lapply(pp, function(x) x$allQps))
qpsByRep <- aggregate(RT~reps*accuracy, qpsByRep, mean)
accuraciesByRep <- do.call(rbind, lapply(pp, function(x) x$accByRep))
accuraciesByRep <- aggregate(accuracy~reps, accuraciesByRep, mean)

# Quantiles in data
dataQps <- lapply(data, getQuantiles, qps=qps)
DqpsCorrect <- do.call(rbind, lapply(dataQps, function(x) x$qpsCorrect))
DqpsError <- do.call(rbind, lapply(dataQps, function(x) x$qpsError))
Daccuracies <- do.call(rbind, lapply(dataQps, function(x) x$accuracy))
points(rep(mean(Daccuracies), 5), apply(DqpsCorrect, 2, mean), pch=4, lwd=2)
points(1-rep(mean(Daccuracies), 5), apply(DqpsError, 2, mean), pch=4, lwd=2)

#### QQ-plot only
# pdf('./model2_QQ.pdf', width=5, height=3, useDingbats = FALSE)
# par(las=1, oma=c(0, 3.5, 0, 0), mar=c(3,0,4,0), mgp=c(2,1,0), xpd=FALSE)
# LO <- cbind(matrix(c(1,1,1,1,1,1,2,1,3,1), nrow=2, ncol=5, byrow=TRUE),
#             matrix(c(1,1,1,1,1,1,2,1,3,1)+3, nrow=2, ncol=5, byrow=TRUE))
# LO <- layout(LO, heights=rep(c(0.1, 3), 1), widths=rep(c(0.2, 3, .1, 3, .2), 2))
# makeQQplotpanel()
# dev.off()

# This is the same data as in variable `data`, but organized differently (useful for plot 2B)
dat$acc <- as.numeric(dat$S)==as.numeric(dat$R)
dat$congruence <- NA
dat$congruence[(dat$cue=="cueleft" & dat$S == 'stimleft') | (dat$cue=="cueright" & dat$S == 'stimright')] <- 'congruent'
dat$congruence[(dat$cue=="cueleft" & dat$S == 'stimright') | (dat$cue=="cueright" & dat$S == 'stimleft')] <- 'incongruent'
dat$congruence[(dat$cue=="cueneutral")] <- 'neutral'
dat$congruence <- factor(dat$congruence, levels=c('congruent', 'neutral', 'incongruent'))

# Start building composite figure -----------------------------------------
# Create layout matrix
tmp <- matrix(1, byrow=TRUE, ncol=100, nrow=100)
tmp[3:47, 52:64] <- 2  ## accuracy
tmp[3:47, 72:84] <- 3
tmp[3:47, 87:99] <- 4
tmp[3:47+50, 53:98-50] <- 5
tmp[3:47+50, 51:75] <- 6

tmp[3:47+50, 52:62] <- 7
tmp[3:47+50, 64:74] <- 8
tmp[3:47+50, 76:100] <- 9
tmp[3:47+50, 77:87] <- 10
tmp[3:47+50, 89:99] <- 11

layout.show(layout(tmp))
layout(tmp)
LO <- cbind(matrix(c(1,1,1,1,1,1,2,1,3,1), nrow=2, ncol=5, byrow=TRUE), 
            matrix(c(1,1,1,1,1,1,2,1,3,1)+3, nrow=2, ncol=5, byrow=TRUE))
layout.show(layout(LO, heights=rep(c(0.1, 3), 1), widths=rep(c(0.2, 3, .1, 3, .2), 2)))

## Get summary statistics of behavioral data: data MRT and acc by difficulty / cue
mrts <- aggregate(RT~coherence*congruence*acc*s, dat, mean)
maccs <- aggregate(acc~coherence*congruence*s, dat, mean)
mrtsM <- aggregate(RT~coherence*congruence*acc, mrts, mean)
maccsM <- aggregate(acc~coherence*congruence, maccs, mean)

# for standard errors, get within-subject variance
mrtBySub <- aggregate(RT~s, dat, mean)
mrtsNormalized <- mrts 
mrtsNormalized$RT <- mrtsNormalized$RT - mrtBySub[as.numeric(as.character(mrts$s)),'RT']
mrtsSE <- aggregate(RT~coherence*congruence*acc, mrts, function(x) sd(x)/sqrt(length(x)))

#
maccBySub <- aggregate(acc~s, dat, mean)
maccsNormalized <- maccs 
maccsNormalized$acc <- maccsNormalized$acc - maccBySub[as.numeric(as.character(maccs$s)),'acc']
maccsSE <- aggregate(acc~coherence*congruence, maccsNormalized, function(x) sd(x)/sqrt(length(x)))

# accuracy
plotLines <- function(means, ses, lty=1) {
  lines(1:3, means, lty=lty)
  points(1:3, means, pch=19)
  arrows(1:3, 
         y0 = means - ses, 
         y1 = means + ses, angle=90, code=3, length=.05) 
  
}


# Composite figure 1 ------------------------------------------------------
#pdf(file='./figure1_big_pt10.pdf', width=6.25, height=4, useDingbats = FALSE, pointsize = 10)
layout(tmp); par(mar=c(2,0,2,0), las=1, mgp=c(2,.75,0), xpd=FALSE)
plot.new()
plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.5, 0.9), main='Accuracy', ylab='Accuracy', xaxt='n')
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
mtext('Accuracy', 2, cex=.66, las=0, line=2)
abline(h=seq(.5, 1, .05), col='lightgrey')
plotLines(means=maccsM[maccsM$coherence=='low','acc'], ses=maccsSE[maccsSE$coherence=='low','acc'])
plotLines(means=maccsM[maccsM$coherence=='high','acc'], ses=maccsSE[maccsSE$coherence=='high','acc'], lty=2)

plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.7, 0.9), main='RT correct', ylab='RT (s)', xaxt='n', yaxt='n')
abline(h=seq(.5, 1, .025), col='lightgrey')
mtext('Mean RT (s)', 2, cex=.66, las=0, line=2)
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
axis(2, at=seq(.7, .9, .05), labels=c('0.7', NA, '0.8', NA, '0.9'))
plotLines(means=mrtsM[mrtsM$coherence=='low'&mrtsM$acc==TRUE,'RT'], ses=mrtsSE[mrtsSE$coherence=='low'&mrtsM$acc==TRUE,'RT'])
plotLines(means=mrtsM[mrtsM$coherence=='high'&mrtsM$acc==TRUE,'RT'], ses=mrtsSE[mrtsSE$coherence=='high'&mrtsM$acc==TRUE,'RT'], lty=2)

plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.7, 0.9), main='RT error', ylab='', xaxt='n', yaxt='n')
abline(h=seq(.5, 1, .025), col='lightgrey')
axis(2, seq(.7, .9, .05), labels=FALSE)
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
plotLines(means=mrtsM[mrtsM$coherence=='low'&mrtsM$acc==FALSE,'RT'], ses=mrtsSE[mrtsSE$coherence=='low'&mrtsM$acc==FALSE,'RT'])
plotLines(means=mrtsM[mrtsM$coherence=='high'&mrtsM$acc==FALSE,'RT'], ses=mrtsSE[mrtsSE$coherence=='high'&mrtsM$acc==FALSE,'RT'], lty=2)
legend('bottomleft', c('Easy', 'Hard'), lty=c(1,2), pch=c(19, 19), bty='n', title='Difficulty')


# BPIC differences --------------------------------------------------------
load('./BPICs.RData')
BPICs <- data.frame(do.call(cbind, lapply(allBPICs, function(x) x[,2])))
colnames(BPICs) <- c('M1', 'M2', 'M3')
sBPICs <- apply(BPICs, 2, sum)
deltaBPICs <- sBPICs-min(sBPICs)
library(gplots)
par(mar=c(2,3,2,2))

# 
deltaBPICs12 <- apply(BPICs, 1, function(x) {x-x[1]})
ordering <- order(deltaBPICs21[2,])
winningModels <- apply(BPICs, 1, which.min)
barplot(deltaBPICs12[2,][ordering], col=winningModels[ordering]+1, xlab='Participant', names=1:length(ordering),
        ylab=expression(paste(Delta, 'BPIC M2 - M1')), main='Model comparison')
legend('bottomright', 
       sapply(c(bquote("M1"~Sigma*"BPIC" == .(round(sBPICs[1],2))),
                bquote("M2"~Sigma*"BPIC" == .(round(sBPICs[2],2))),
                bquote("M3"~Sigma*"BPIC" == .(round(sBPICs[3],2)))), as.expression),
       bty='n', col=2:4, pch=15, title = 'Winning model')


par(mar=c(2,0,2,0))
makeQQplotpanel()
dev.off()



# Composite figure 2, with Bayesian Averaged Parameter estimates --------------------------------------------
bmaPars <- read.csv('../derivatives/BMAParameters.csv')[,2:9]
bmaPars <- bmaPars[,c('vshift', 'vshiftBias', 'zshiftBias')]
means <- apply(bmaPars, 2, mean)
se <- apply(bmaPars, 2, function(x) sd(x)/sqrt(length(x)))

# par(las=1)
# barplot2(means, ci.l = means-se, ci.u=means+se, plot.ci=TRUE, xlab='Parameter value',
#          names=c(expression(paste(Delta, 'v'['difficulty'])),
#                  expression(paste(Delta, 'v'['cue'])),
#                  expression(paste(Delta, 'z'['cue']))),
#          horiz=FALSE)


# Create new layout --------------------------------------------------------
tmp <- matrix(1, byrow=TRUE, ncol=100, nrow=100)
tmp[3:47, 52:64] <- 2  ## accuracy
tmp[3:47, 72:84] <- 3
tmp[3:47, 87:99] <- 4
tmp[3:47+50, 53:82-50] <- 5
tmp[3:47+50, 83:98-50] <- 6
tmp[3:47+50, 51:75] <- 7
tmp[3:47+50, 52:62] <- 8
tmp[3:47+50, 64:74] <- 9
tmp[3:47+50, 76:100] <- 10
tmp[3:47+50, 77:87] <- 11
tmp[3:47+50, 89:99] <- 12
layout.show(layout(tmp))

# Start plotting
pdf(file='./figure1_big_pt10-2.pdf', width=6.25, height=4, useDingbats = FALSE, pointsize = 10)
layout(tmp); par(mar=c(2,0,2,0), las=1, mgp=c(2,.75,0), xpd=FALSE)
plot.new()
plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.5, 0.9), main='Accuracy', ylab='Accuracy', xaxt='n')
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
mtext('Accuracy', 2, cex=.66, las=0, line=2)
abline(h=seq(.5, 1, .05), col='lightgrey')
plotLines(means=maccsM[maccsM$coherence=='low','acc'], ses=maccsSE[maccsSE$coherence=='low','acc'])
plotLines(means=maccsM[maccsM$coherence=='high','acc'], ses=maccsSE[maccsSE$coherence=='high','acc'], lty=2)
#legend('topright', c('Easy', 'Hard'), lty=c(1,2), pch=c(19, 19), bty='n', title='Difficulty')

plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.7, 0.9), main='RT correct', ylab='RT (s)', xaxt='n', yaxt='n')
abline(h=seq(.5, 1, .025), col='lightgrey')
mtext('Mean RT (s)', 2, cex=.66, las=0, line=2)
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
axis(2, at=seq(.7, .9, .05), labels=c('0.7', NA, '0.8', NA, '0.9'))
plotLines(means=mrtsM[mrtsM$coherence=='low'&mrtsM$acc==TRUE,'RT'], ses=mrtsSE[mrtsSE$coherence=='low'&mrtsM$acc==TRUE,'RT'])
plotLines(means=mrtsM[mrtsM$coherence=='high'&mrtsM$acc==TRUE,'RT'], ses=mrtsSE[mrtsSE$coherence=='high'&mrtsM$acc==TRUE,'RT'], lty=2)
#legend('topright', c('Easy', 'Hard'), lty=c(1,2), pch=c(19, 19), bty='n', title='Difficulty')

plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.7, 0.9), main='RT error', ylab='', xaxt='n', yaxt='n')
abline(h=seq(.5, 1, .025), col='lightgrey')
axis(2, seq(.7, .9, .05), labels=FALSE)
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
plotLines(means=mrtsM[mrtsM$coherence=='low'&mrtsM$acc==FALSE,'RT'], ses=mrtsSE[mrtsSE$coherence=='low'&mrtsM$acc==FALSE,'RT'])
plotLines(means=mrtsM[mrtsM$coherence=='high'&mrtsM$acc==FALSE,'RT'], ses=mrtsSE[mrtsSE$coherence=='high'&mrtsM$acc==FALSE,'RT'], lty=2)
legend('bottomleft', c('Easy', 'Hard'), lty=c(1,2), pch=c(19, 19), bty='n', title='Difficulty')


# BPIC differences --------------------------------------------------------
load('./BPICs.RData')
BPICs <- data.frame(do.call(cbind, lapply(allBPICs, function(x) x[,2])))
colnames(BPICs) <- c('M1', 'M2', 'M3')
sBPICs <- apply(BPICs, 2, sum)
deltaBPICs <- sBPICs-min(sBPICs)
library(gplots)
par(mar=c(2,3,2,1))

# 
deltaBPICs12 <- apply(BPICs, 1, function(x) {x-x[1]})
ordering <- order(deltaBPICs12[2,])
winningModels <- apply(BPICs, 1, which.min)
barplot(deltaBPICs12[2,][ordering], col=winningModels[ordering]+1, xlab='Participant', names=1:length(ordering),
        ylab=expression(paste(Delta, 'BPIC M2 - M1')), main='Model comparison', names.arg = c(''))
legend('bottomright', 
       sapply(c(bquote("M1"~Sigma*"BPIC" == .(round(sBPICs[1]))),
                bquote("M2"~Sigma*"BPIC" == .(round(sBPICs[2]))),
                bquote("M3"~Sigma*"BPIC" == .(round(sBPICs[3])))), as.expression),
       bty='n', col=2:4, pch=15, title = 'Winning model')

par(mar=c(3,3,2,2))
barplot2(means, ci.l = means-se, ci.u=means+se, plot.ci=TRUE, ylab='Parameter value', main='Parameters',
         names=c(expression(paste(Delta, 'v'['difficulty'])),
                 expression(paste(Delta, 'v'['cue'])),
                 expression(paste(Delta, 'z'['cue']))),
         horiz=FALSE, las=2)

par(mar=c(2,0,2,0))
makeQQplotpanel()
dev.off()



#  Composite figure 3, same as above, but switch DDM names ----------------
tmp <- matrix(1, byrow=TRUE, ncol=100, nrow=100)
tmp[3:47, 53:98-50] <- 2
tmp[3:47, 52:64] <- 3  ## accuracy
tmp[3:47, 72:84] <- 4
tmp[3:47, 87:99] <- 5
tmp[3:47+50, 53:82-50] <- 6
tmp[3:47+50, 83:98-50] <- 7
tmp[3:47+50, 51:75] <- 8
tmp[3:47+50, 52:62] <- 9
tmp[3:47+50, 64:74] <- 10
tmp[3:47+50, 76:100] <- 11
tmp[3:47+50, 77:87] <- 12
tmp[3:47+50, 89:99] <- 13
layout.show(layout(tmp))

# Start plottin
## Plot
pdf(file='./figure1_big_pt10-4.pdf', width=6.25, height=4, useDingbats = FALSE, pointsize = 10)
layout(tmp); par(mar=c(2,0,2,0), las=1, mgp=c(2,.75,0), xpd=FALSE)
plot.new(); plot.new(); fig_label('A', cex=1.2, font=2) #text(0, 1, 'A', cex=1)
plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.5, 0.9), main='Accuracy', ylab='Accuracy', xaxt='n')
fig_label('B', cex=1.2, font=2)
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
mtext('Accuracy', 2, cex=.66, las=0, line=2)
abline(h=seq(.5, 1, .05), col='lightgrey')
plotLines(means=maccsM[maccsM$coherence=='low','acc'], ses=maccsSE[maccsSE$coherence=='low','acc'])
plotLines(means=maccsM[maccsM$coherence=='high','acc'], ses=maccsSE[maccsSE$coherence=='high','acc'], lty=2)
#legend('topright', c('Easy', 'Hard'), lty=c(1,2), pch=c(19, 19), bty='n', title='Difficulty')

plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.7, 0.9), main='RT correct', ylab='RT (s)', xaxt='n', yaxt='n')
abline(h=seq(.5, 1, .025), col='lightgrey')
mtext('Mean RT (s)', 2, cex=.66, las=0, line=2)
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
axis(2, at=seq(.7, .9, .05), labels=c('0.7', NA, '0.8', NA, '0.9'))
plotLines(means=mrtsM[mrtsM$coherence=='low'&mrtsM$acc==TRUE,'RT'], ses=mrtsSE[mrtsSE$coherence=='low'&mrtsM$acc==TRUE,'RT'])
plotLines(means=mrtsM[mrtsM$coherence=='high'&mrtsM$acc==TRUE,'RT'], ses=mrtsSE[mrtsSE$coherence=='high'&mrtsM$acc==TRUE,'RT'], lty=2)
#legend('topright', c('Easy', 'Hard'), lty=c(1,2), pch=c(19, 19), bty='n', title='Difficulty')

plot(0,0, type='n', xlim=c(.8, 3.2), xlab='Cue', ylim=c(0.7, 0.9), main='RT error', ylab='', xaxt='n', yaxt='n')
abline(h=seq(.5, 1, .025), col='lightgrey')
axis(2, seq(.7, .9, .05), labels=FALSE)
axis(1, 1:3, labels=c('C', 'N', 'I'), mgp=c(2,.75,0))
mtext('Cue', side=1, line=2, cex=par()$cex)
plotLines(means=mrtsM[mrtsM$coherence=='low'&mrtsM$acc==FALSE,'RT'], ses=mrtsSE[mrtsSE$coherence=='low'&mrtsM$acc==FALSE,'RT'])
plotLines(means=mrtsM[mrtsM$coherence=='high'&mrtsM$acc==FALSE,'RT'], ses=mrtsSE[mrtsSE$coherence=='high'&mrtsM$acc==FALSE,'RT'], lty=2)
legend('bottomleft', c('Easy', 'Hard'), lty=c(1,2), pch=c(19, 19), bty='n', title='Difficulty')


# BPIC differences --------------------------------------------------------
load('./BPICs.RData')
BPICs <- data.frame(do.call(cbind, lapply(allBPICs, function(x) x[,2])))
colnames(BPICs) <- c('M1', 'M2', 'M3')
sBPICs <- apply(BPICs, 2, sum)
deltaBPICs <- sBPICs-min(sBPICs)
library(gplots)
par(mar=c(2,3,2,1))

# 
deltaBPICs12 <- apply(BPICs, 1, function(x) {x-x[1]})
ordering <- order(deltaBPICs12[2,])
winningModels <- apply(BPICs, 1, which.min)
barplot(deltaBPICs12[2,][ordering], col=winningModels[ordering]+1, xlab='Participant', names=1:length(ordering),
        ylab=expression(paste(Delta, 'BPIC M2 - M1')), main='Model comparison', names.arg = c(''))
legend('bottomright', 
       sapply(c(bquote("M1"~Sigma*"BPIC" == .(round(sBPICs[1]))),
                bquote("M2"~Sigma*"BPIC" == .(round(sBPICs[2]))),
                bquote("M3"~Sigma*"BPIC" == .(round(sBPICs[3])))), as.expression),
       bty='n', col=2:4, pch=15, title = 'Winning model')
fig_label('C', cex=1.2, font=2)

par(mar=c(3,3,2,2))
barplot2(means, ci.l = means-se, ci.u=means+se, plot.ci=TRUE, ylab='Parameter value', main='Parameters',
         names=c(expression(paste(Delta, 'v'['difficulty'])),
                 expression(paste(Delta, 'v'['cue'])),
                 expression(paste(Delta, 'z'['cue']))),
         horiz=FALSE, las=2)
fig_label('D', cex=1.2, font=2)
par(mar=c(2,0,2,0))
makeQQplotpanel('E')
dev.off()

