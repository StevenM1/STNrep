allBPICs <- list()
allSummaries <- list()
for(modelN in c(1,2,3)) {
  source(paste0("model", modelN, "sv.R"))
  data.model <- data.model.dmc(dat, model)
  load(paste0('samplesModel', modelN, 'sv.RData'))
  #  samples <- do.load.samples(modelN, data.model, p.prior)
  
  allBPICs[[modelN]] = h.IC.dmc(samples)
  allSummaries[[modelN]] = summary.dmc(samples)
}

# Model comparisons -------------------------------------------------------
# 1. Sum BPICs across participants
lapply(allBPICs, function(x) apply(x, 2, sum))
# model 3 wins if we sum all BPICs measures, however....

BPICbySub <- do.call(cbind, lapply(allBPICs, function(x) x[,2]))
#allBPICs <- cbind(BPIC1[,2], BPIC2[,2], BPIC3[,2], BPIC4[,2])
#allBPICs <- cbind(BPIC2[,2], BPIC3[,2], BPIC4[,2])
table(apply(BPICbySub, 1, which.min))
plot(table(apply(BPICbySub, 1, which.min)))
# ...model 1 wins for most participants individually. Bayesian model averaging? Get wBPIC

minBPICs = apply(BPICbySub, 1, min)
deltaBPICs = BPICbySub - minBPICs
wBPICs = exp(-.5 * (deltaBPICs)) / apply(exp(-.5 * (deltaBPICs)), 1, sum)
apply(wBPICs, 2, mean)  
# model 3 wins again by wBPIC - it appears model 1 fits best for most participants, but fits really poorly for some others

medianParams <- lapply(allSummaries, function(x) data.frame(do.call(rbind, lapply(x, function(y) y$quantiles[,3]))))
medianParams[[1]]$vshift = medianParams[[1]]$v.high - medianParams[[1]]$v.low
medianParams[[2]]$vshift = medianParams[[2]]$v.high - medianParams[[2]]$v.low
medianParams[[3]]$vshift = medianParams[[3]]$v.high - medianParams[[3]]$v.low

# For model 1, add drift shift of bias (set to 0)
medianParams[[1]]$vshiftBias = 0
medianParams[[2]]$vshiftBias = medianParams[[2]]$v_bias.biastrial
medianParams[[3]]$vshiftBias = medianParams[[3]]$v_bias.biastrial

# for model 2, add start point shift of bias (set to 0)
medianParams[[1]]$zshiftBias = medianParams[[1]]$z_bias.biastrial
medianParams[[2]]$zshiftBias = 0
medianParams[[3]]$zshiftBias = medianParams[[3]]$z_bias.biastrial

# BMA weighting
parameters = c('a', 'v.low', 'v.high', 'z', 'zshiftBias', 'vshift', 'vshiftBias', 't0')
bmaParameters = medianParams[[1]][,parameters] * wBPICs[,1] + 
  medianParams[[2]][,parameters] * wBPICs[,2] + 
  medianParams[[3]][,parameters] * wBPICs[,3]


bmaParameters.withsv <- bmaParameters
bmaParameters.nosv <- read.csv('../derivatives/BMAParameters.csv')

##
par(mfrow=c(3,3))
for(parName in colnames(bmaParameters.withsv)) {
  plot(bmaParameters.withsv[,parName], bmaParameters.nosv[,parName], xlab='With sv and sz', ylab='Without sv and sz', main=parName)
  legend('bottomright', c(paste0('r = ', round(cor(bmaParameters.withsv[,parName], bmaParameters.nosv[,parName]), 2))), lty=c(1), pch=c(NA), col=c(NA), bty='n')
}

## almost exactly the same vdifference, zshift and vshit