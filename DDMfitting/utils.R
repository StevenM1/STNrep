runDiagnostics <- function(fn, dat, samples=NULL, nCores=30, verbose=TRUE, chainsFileType='jpeg',
                           plotChains=TRUE, plotPosteriorPriors=TRUE, calculateESS=TRUE, plotFits=TRUE, calculateParametersBPIC=TRUE) {
  if(is.null(samples)) samples <- loadSamples(fn, samplesDir='samples')
  
  saveDir <- file.path('diagnostics', fn)
  dir.create(saveDir, showWarnings = FALSE)
  if(chainsFileType == 'jpeg') dir.create(file.path(saveDir, 'chains'), showWarnings = FALSE)
  
  # Chains: Plot hyper, append Gelman's
  if(plotChains) {
    if(verbose) print('Plotting chains...')
    gms <- h.gelman.diag.dmc(samples)
    if(chainsFileType == 'pdf') {
      pdf(file=file.path(saveDir, 'chains.pdf'))
    } else {
      jpeg(filename=file.path(saveDir, 'chains', 'chains-%03d.jpeg'), width=7, height=7, units='in', quality=100, res=200)
    }
    par(mpg=c(1,1,0), mar=c(4,3,3,1)+.1, oma=c(0,0,1,0), mfrow=c(4,4))
#    plot.dmc(samples, hyper=TRUE, layout=c(4,4))
#    title(paste0('Hyper, Gelmans diag = ', round(gms['hyper'], 4)), outer=TRUE, line=0)
    
    # Chains: Plot per participant, append Gelman's
    for(subject in names(samples)) {
      plot.dmc(samples, subject=subject, layout=c(4,4)); 
      title(paste0('Subject ', subject, ' Gelmans diag = ', round(gms[which(names(gms)==as.character(subject))], 4)), outer=TRUE, line=0)
    }
    dev.off()
  }
  
  # Posterior vs Prior
  if(plotPosteriorPriors) {
    if(!grepl('softmax', fn)) {
      pdf(file=file.path(saveDir, 'posteriorPrior.pdf'))
      if(verbose) print('Plotting posterior vs prior...')
      par(mpg=c(1,1,0), mar=c(4,3,3,1)+.1, oma=c(0,0,1,0), mfrow=c(4,4))
#      plot.dmc(samples, hyper=TRUE, p.prior = attr(samples, 'hyper')$pp.prior, layout=c(4,4))
      for(subject in names(samples)) {
        plot.dmc(samples, subject=subject, p.prior = samples[[subject]]$p.prior, layout=c(4,4))
        title(paste0('Subject ', subject, ' Gelmans diag = ', round(gms[which(names(gms)==as.character(subject))], 4)), outer=TRUE, line=0)
      } 
      dev.off()
    } else {
      warning('Skipping posterior vs prior plot, this crashes for softmax... (too few parameters?)')
    }
  }
  
  # Effective sample size: save to csv
  if(calculateESS) {
    if(verbose) print('Calculating effective samples sizes...')
    effS <- data.frame(do.call(rbind, effectiveSize.dmc(samples)))
#    effSH <- effectiveSize.dmc(samples, hyper=TRUE)
#    effS <- rbind(effS, effSH[grepl('.h1', names(effSH))])
#    effS <- rbind(effS, effSH[grepl('.h2', names(effSH))])
    effS <- rbind(effS, apply(effS, 2, min))
    row.names(effS) <- c(row.names(effS)[1:(nrow(effS)-1)], 'minimum')
    write.csv(effS, file=file.path(saveDir, 'effectiveSizes.csv'))
  }
  
  # Fits
  if(plotFits) {
    if(verbose) print('Plotting posterior predictives...')
    pdf(file=file.path(saveDir, 'posteriorPredictives.pdf'))
    data <- lapply(samples, function(x) x$data)
    ppSim <- h.post.predict.dmc(samples, cores=30, save.simulation = TRUE)
    ppSim <- lapply(ppSim, getAccuracy)
    data <- lapply(data, getAccuracy)
    do.plot.overall(ppSim, data, qps=seq(.01, .99, .01), model=gsub('model', '', x = fn))
    dev.off()
  }
  
  # Parameters: save to csv, append BPIC per participant
  if(calculateParametersBPIC) {
    if(verbose) print('Getting parameters, BPIC...')
    bpics <- h.IC.dmc(samples)
    summ <- summary.dmc(samples)
    medians <- data.frame(do.call(rbind, lapply(summ, function(x) x$quantiles[,3])))
    medians <- rbind(medians, apply(medians, 2, mean))
    medians <- rbind(medians, apply(medians, 2, sd))
    
    medians$minimum.deviances <- c(bpics[,1], sum(bpics[,1]), NA)
    medians$BPIC <- c(bpics[,2], sum(bpics[,2]), NA)
    row.names(medians) <- c(row.names(medians)[1:(nrow(medians)-2)], 'mean/sum', 'sd')
    write.csv(medians, file=file.path(saveDir, 'medianParametersBPICs.csv'))
  }
}




## 'better' plotting: Ignore left/right stim, plot accuracy-coded dCDFs
do.plot.quantiles <- function(pp, data, qps=seq(.01, .99, .01), main, xlim=c(0.4, 2.5), ylim=c(0, .9)) {
  pp <- lapply(pp, getQuantiles, qps=qps)
  
  qpsCorrect <- do.call(rbind, lapply(pp, function(x) x$qpsCorrect))
  qpsError <- do.call(rbind, lapply(pp, function(x) x$qpsError))
  accuracies <- do.call(rbind, lapply(pp, function(x) x$accuracy))
  qpsByRep <- do.call(rbind, lapply(pp, function(x) x$allQps))
  qpsByRep <- aggregate(RT~reps*accuracy, qpsByRep, mean)
  accuraciesByRep <- do.call(rbind, lapply(pp, function(x) x$accByRep))
  accuraciesByRep <- aggregate(accuracy~reps, accuraciesByRep, mean)
  
  ## Model
  plot(0, 0, xlim=xlim, ylim=ylim, type='n', xlab='RT', ylab='Def CDF', main=main)
  for(rep in qpsByRep$reps) {
    points(qpsByRep[qpsByRep$accuracy==1&qpsByRep$reps==rep,seq(10, 90, 20)+2],
           qps[seq(10, 90, 20)]*accuraciesByRep[accuraciesByRep$reps==rep,'accuracy'], pch=20, col='grey')
    points(qpsByRep[qpsByRep$accuracy==0&qpsByRep$reps==rep,seq(10, 90, 20)+2],
           qps[seq(10, 90, 20)]*(1-accuraciesByRep[accuraciesByRep$reps==rep,'accuracy']), pch=20, col='grey')
  }
  lines(apply(qpsCorrect, 2, mean), qps*mean(accuracies), lwd=1)
  lines(apply(qpsError, 2, mean), qps*(1-mean(accuracies)), lwd=1)
  points(apply(qpsCorrect, 2, mean)[seq(10, 90, 20)], qps[seq(10, 90, 20)]*mean(accuracies))
  points(apply(qpsError, 2, mean)[seq(10, 90, 20)], qps[seq(10, 90, 20)]*(1-mean(accuracies)))
  
  ## Data
  dataQps <- lapply(data, getQuantiles, qps=qps)
  DqpsCorrect <- do.call(rbind, lapply(dataQps, function(x) x$qpsCorrect))
  DqpsError <- do.call(rbind, lapply(dataQps, function(x) x$qpsError))
  Daccuracies <- do.call(rbind, lapply(dataQps, function(x) x$accuracy))
  
  lines(apply(DqpsCorrect, 2, mean), qps*mean(Daccuracies), lty=2, col='red', lwd=2)
  lines(apply(DqpsError, 2, mean), qps*(1-mean(Daccuracies)), lty=2, col='red', lwd=2)
  points(apply(DqpsCorrect, 2, mean)[seq(10, 90, 20)], qps[seq(10, 90, 20)]*mean(Daccuracies), col='red', lwd=2)
  points(apply(DqpsError, 2, mean)[seq(10, 90, 20)], qps[seq(10, 90, 20)]*(1-mean(Daccuracies)), col='red', lwd=2)
  legend('topright', c('Data', 'Model'), col=c('red', 'black'), lty=c(2, 1), bty='n', lwd=c(2,1))
}

do.plot.overall <- function(pp, data, model, qps, ...) {
  ## plot overall
  par(mfrow=c(2,3))
  do.plot.quantiles(pp, data=data, qps=qps, main=paste0('Model ', model, ', overall'))
  
  ## Plot by difficulty
  ppEasy <- lapply(pp, function(x) x[x$coherence=='high',])
  dataEasy <- lapply(data, function(x) x[x$coherence=='high',])
  do.plot.quantiles(ppEasy, dataEasy, qps=qps, main=paste0('Model ', model, ', easy'))
  
  ppHard <- lapply(pp, function(x) x[x$coherence=='low',])
  dataHard <- lapply(data, function(x) x[x$coherence=='low',])
  do.plot.quantiles(ppHard, dataHard, qps=qps, main=paste0('Model ', model, ', hard'))
  
  ## Plot by cue congruency
  ppCongruent <- lapply(pp, function(x) x[x$S=='stimleft' & x$cue=='cueleft' | x$S=='stimright' & x$cue=='cueright',])
  dataCongruent <- lapply(data, function(x) x[x$S=='stimleft' & x$cue=='cueleft' | x$S=='stimright' & x$cue=='cueright',])
  ppIncongruent <- lapply(pp, function(x) x[x$S=='stimleft' & x$cue=='cueright' | x$S=='stimright' & x$cue=='cueleft',])
  dataIncongruent <- lapply(data, function(x) x[x$S=='stimleft' & x$cue=='cueright' | x$S=='stimright' & x$cue=='cueleft',])
  ppNeutral <- lapply(pp, function(x) x[x$cue=='cueneutral',])
  dataNeutral <- lapply(data, function(x) x[x$cue=='cueneutral',])
  
  do.plot.quantiles(ppCongruent, dataCongruent, qps=qps, main=paste0('Model ', model, ', congruent'))
  do.plot.quantiles(ppIncongruent, dataIncongruent, qps=qps, main=paste0('Model ', model, ', incongruent'))
  do.plot.quantiles(ppNeutral, dataNeutral, qps=qps, main=paste0('Model ', model, ', neutral'))
}

getAccuracy <- function(x) {
  x$accuracy <- as.numeric(x$S) == as.numeric(x$R)
  return(x)
}
getQuantiles <- function(x, qps) {
  if(!'reps' %in% colnames(x)) x$reps <- 1
  allQps <- aggregate(RT~accuracy*reps, x, quantile, probs=qps)
  meanQps <- aggregate(RT~accuracy, allQps, mean)
  
  #  qpsCorrect <- quantile(x$RT[x$accuracy==1], probs=qps)
  #  qpsError <- quantile(x$RT[x$accuracy==0], probs=qps)
  accByRep <- aggregate(accuracy~reps, x, mean)
  meanAccuracy <- mean(x$accuracy)
  return(list('qpsCorrect'=meanQps[meanQps$accuracy==1,2:(length(qps)+1)], 
              'qpsError'=meanQps[meanQps$accuracy==0,2:(length(qps)+1)], 
              'accuracy'=meanAccuracy,
              'allQps'=allQps,
              'accByRep'=accByRep))
}
