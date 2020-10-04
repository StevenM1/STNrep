do.load.samples <- function(modelN, data.model, p.prior, cores=30) {
  samplesFn <- paste0('./samplesModel', modelN, '.RData')
  if(file.exists(samplesFn)) {
    load(samplesFn)
#    samples <- h.samples.dmc(nmc=0, add=TRUE, samples=samples, remove=1:500)
  } else {
    print('Starting burn-in...')
    # simple burn-in: 500 samples will suffice
    samples <- h.samples.dmc(nmc = 500, p.prior, data.model, thin = 1)
    samples <- h.run.dmc(samples, cores=cores, p.migrate = 0.05)
    
    print('Starting sampling from posterior...')
    # add 5000 samples at a time, using DMCs h.RUN.dmc
    samples <- h.samples.dmc(nmc = 5000, samples=samples, add=FALSE)
    samples <- h.RUN.dmc(samples, cores=cores)
    save(samples, file=samplesFn)
  }
  return(samples)
}

do.make.traces <- function(modelN, samples, check.exists=TRUE) {
  tracesFn <- paste0('./tracesModel', modelN, '.pdf')
  if(!check.exists | (check.exists & !file.exists(tracesFn))) {
    pdf(tracesFn, width=9, height=9)
    for(i in 1:length(samples)) {
      plot.dmc(samples[[i]], layout=c(3,3))
      mtext(i, side=1, outer=TRUE, line=0)
    }
    dev.off()
  }
}


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
do.plot.overall <- function(pp, data, model, qps, subject=NULL, ...) {
  if(!is.null(subject)) {
    pp <- list(pp[[subject]])
    data <- list(data[[subject]])
    titleText <- paste0('Model ', modelN, ', subject ', subject)
  } else {
    titleText <- paste0('Model ', modelN, ', across subjects')
  }
  ## plot overall
  par(mfrow=c(2,3), oma=c(1,0,1,1)+.1)
  do.plot.quantiles(pp, data=data, qps=qps, main='Overall')
  
  ## Plot by difficulty
  ppEasy <- lapply(pp, function(x) x[x$coherence=='high',])
  dataEasy <- lapply(data, function(x) x[x$coherence=='high',])
  do.plot.quantiles(ppEasy, dataEasy, qps=qps, main='Easy')
  
  ppHard <- lapply(pp, function(x) x[x$coherence=='low',])
  dataHard <- lapply(data, function(x) x[x$coherence=='low',])
  do.plot.quantiles(ppHard, dataHard, qps=qps, main='Hard')
  
  ## Plot by cue congruency
  ppCongruent <- lapply(pp, function(x) x[x$S=='stimleft' & x$cue=='cueleft' | x$S=='stimright' & x$cue=='cueright',])
  dataCongruent <- lapply(data, function(x) x[x$S=='stimleft' & x$cue=='cueleft' | x$S=='stimright' & x$cue=='cueright',])
  ppIncongruent <- lapply(pp, function(x) x[x$S=='stimleft' & x$cue=='cueright' | x$S=='stimright' & x$cue=='cueleft',])
  dataIncongruent <- lapply(data, function(x) x[x$S=='stimleft' & x$cue=='cueright' | x$S=='stimright' & x$cue=='cueleft',])
  ppNeutral <- lapply(pp, function(x) x[x$cue=='cueneutral',])
  dataNeutral <- lapply(data, function(x) x[x$cue=='cueneutral',])
  
  do.plot.quantiles(ppCongruent, dataCongruent, qps=qps, main='Congruent')
  do.plot.quantiles(ppIncongruent, dataIncongruent, qps=qps, main='Incongruent')
  do.plot.quantiles(ppNeutral, dataNeutral, qps=qps, main='Neutral')
  mtext(titleText, side=3, line=-1, outer=TRUE)
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


getPlottingValues <- function(pp, data) {
  pp <- lapply(pp, getQuantiles, qps=qps)
  
  qpsCorrect <- do.call(rbind, lapply(pp, function(x) x$qpsCorrect))
  qpsError <- do.call(rbind, lapply(pp, function(x) x$qpsError))
  accuracies <- do.call(rbind, lapply(pp, function(x) x$accuracy))
  qpsByRep <- do.call(rbind, lapply(pp, function(x) x$allQps))
  qpsByRep <- aggregate(RT~reps*accuracy, qpsByRep, mean)
  accuraciesByRep <- do.call(rbind, lapply(pp, function(x) x$accByRep))
  accuraciesByRep <- aggregate(accuracy~reps, accuraciesByRep, mean)
  
  dataQps <- lapply(data, getQuantiles, qps=qps)
  DqpsCorrect <- do.call(rbind, lapply(dataQps, function(x) x$qpsCorrect))
  DqpsError <- do.call(rbind, lapply(dataQps, function(x) x$qpsError))
  Daccuracies <- do.call(rbind, lapply(dataQps, function(x) x$accuracy))
  
  return(list('qpsCorrect'=qpsCorrect,
              'qpsError'=qpsError,
              'accuracies'=accuracies,
              'qpsByRep'=qpsByRep,
              'accuraciesByRep'=accuraciesByRep,
              'DqpsCorrect'=DqpsCorrect,
              'DqpsError'=DqpsError,
              'Daccuracies'=Daccuracies))
}

plotQQ <- function(plottingValues, colData=1, colModel=1) {
  list2env(plottingValues, globalenv())
  
  for(rep in unique(qpsByRep$reps)) {
    thisAcc <- accuraciesByRep[rep,'accuracy']
    points(rep(thisAcc, 5), qpsByRep[qpsByRep$accuracy==1&qpsByRep$reps==rep,3:7], pch=20, col='grey')
    points(rep(1-thisAcc, 5), qpsByRep[qpsByRep$accuracy==0&qpsByRep$reps==rep,3:7], pch=20, col='grey')
  }
  # points(rep(mean(accuracies), 5), apply(qpsCorrect, 2, mean), lwd=2, pch=16, col=colModel)
  # points(1-rep(mean(accuracies), 5), apply(qpsError, 2, mean), lwd=2, pch=16, col=colModel)
  # 
  # points(rep(mean(Daccuracies), 5), apply(DqpsCorrect, 2, mean), pch=4, lwd=2, col=colData)
  # points(1-rep(mean(Daccuracies), 5), apply(DqpsError, 2, mean), pch=4, lwd=2, col=colData)
  
  points(rep(mean(accuracies), 5), apply(qpsCorrect, 2, mean), lwd=2.1, pch=16, col=1)
  points(1-rep(mean(accuracies), 5), apply(qpsError, 2, mean), lwd=2.1, pch=16, col=1)
  
  points(rep(mean(Daccuracies), 5), apply(DqpsCorrect, 2, mean), pch=4, lwd=2.1, col=1)
  points(1-rep(mean(Daccuracies), 5), apply(DqpsError, 2, mean), pch=4, lwd=2.1, col=1)
  
  points(rep(mean(accuracies), 5), apply(qpsCorrect, 2, mean), lwd=2, pch=16, col=colModel)
  points(1-rep(mean(accuracies), 5), apply(qpsError, 2, mean), lwd=2, pch=16, col=colModel)
  
  points(rep(mean(Daccuracies), 5), apply(DqpsCorrect, 2, mean), pch=4, lwd=2, col=colData)
  points(1-rep(mean(Daccuracies), 5), apply(DqpsError, 2, mean), pch=4, lwd=2, col=colData)
}

makeQQplotpanel <- function(addFigLabel=NULL) {
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  if(!is.null(addFigLabel)) fig_label(addFigLabel, cex=1.2, font=2)
  mtext('By difficulty', side=3, cex=par()$cex.main*par()$cex, font=2, line=1)
  mtext('Choice proportion', side=1, cex=par()$cex, line=2)
  
  plot(0, 0, type='n', xlim=c(0.15, .4), ylim=c(0.5, 1.5), ylab='', xlab='', xaxt='n')
  mtext('Error', line=0, cex=par()$cex, font=2)
  axis(1, at=seq(.2, .4, .1))
  mtext('RT (s)', side=2, line=2, cex=.66, las=0)
  abline(h=seq(0, 2, .1), col='lightgrey')
  abline(v=seq(0, 2, .1), col='lightgrey')
  ppEasy <- lapply(ppSim, function(x) x[x$coherence=='high',])
  dataEasy <- lapply(data, function(x) x[x$coherence=='high',])
  plotQQ(getPlottingValues(ppEasy, dataEasy))
  legend('topleft', c('D', 'M'), col=c(1,1), pch=c(4,16), lwd=2, bty='n', x.intersp = 0, lty=NA)
  
  ppHard <- lapply(ppSim, function(x) x[x$coherence=='low',])
  dataHard <- lapply(data, function(x) x[x$coherence=='low',])
  plotQQ(getPlottingValues(ppHard, dataHard), colData=2, colModel=2)
  
  #
  plot(0, 0, type='n', xlim=c(0.65, .85), ylim=c(0.5, 1.5), main='', ylab='', xlab='', yaxt='n', xaxt='n;')
  mtext('Correct', line=0, cex=par()$cex, font=2)
  axis(1, at=seq(.7, .8, .1))
  #axis(2, at=seq(.6, 1.4, .2), labels=FALSE)
  abline(h=seq(0, 2, .1), col='lightgrey')
  abline(v=seq(0, 2, .1), col='lightgrey')
  plotQQ(getPlottingValues(ppEasy, dataEasy))
  plotQQ(getPlottingValues(ppHard, dataHard), colData=2, colModel=2)
  legend('topright', c('E', 'H'), col=c(1,2), pch=c(15,15), lwd=2, bty='n', x.intersp = 0, lty=NA)
  
  
  ##### cues
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  mtext('By cue', side=3, cex=par()$cex.main*par()$cex, font=2, line=1)
  mtext('Choice proportion', side=1, cex=par()$cex, line=2)
  
  plot(0, 0, type='n', xlim=c(0.15, .4), ylim=c(0.5, 1.5), main='', ylab='', xlab='', yaxt='n', xaxt='n')
  mtext('Error', line=0, cex=par()$cex, font=2)
  axis(1, at=seq(.2, .4, .1))
  axis(2, at=seq(.6, 1.4, .2), labels=FALSE)
  abline(h=seq(0, 2, .1), col='lightgrey')
  abline(v=seq(0, 2, .1), col='lightgrey')
  ppCongruent <- lapply(ppSim, function(x) x[x$S=='stimleft' & x$cue=='cueleft' | x$S=='stimright' & x$cue=='cueright',])
  ppIncongruent <- lapply(ppSim, function(x) x[x$S=='stimleft' & x$cue=='cueright' | x$S=='stimright' & x$cue=='cueleft',])
  ppNeutral <- lapply(ppSim, function(x) x[x$cue=='cueneutral',])
  
  dataIncongruent <- lapply(data, function(x) x[x$S=='stimleft' & x$cue=='cueright' | x$S=='stimright' & x$cue=='cueleft',])
  dataCongruent <- lapply(data, function(x) x[x$S=='stimleft' & x$cue=='cueleft' | x$S=='stimright' & x$cue=='cueright',])
  dataNeutral <- lapply(data, function(x) x[x$cue=='cueneutral',])
  
  plotQQ(getPlottingValues(ppCongruent, dataCongruent))
  plotQQ(getPlottingValues(ppIncongruent, dataIncongruent), colModel=3, colData=3)
  plotQQ(getPlottingValues(ppNeutral, dataNeutral), colModel=2, colData=2)
  legend('topleft', c('D', 'M'), col=c(1,1), pch=c(4,16), lwd=2, bty='n', x.intersp = 0, lty=NA)
  
  #par(las=1, oma=c(0, 2, 0, 0), mar=c(0,2,3,0), mgp=c(2,1,0), xpd=FALSE)
  plot(0, 0, type='n', xlim=c(0.6, .85), ylim=c(0.5, 1.5), main='', ylab='', xlab='', yaxt='n', xaxt='n')
  axis(1, at=seq(.7, .8, .1))
  mtext('Correct', line=0, cex=par()$cex, font=2)
  #axis(2, at=seq(.6, 1.4, .2), labels=FALSE)
  abline(h=seq(0, 2, .1), col='lightgrey')
  abline(v=seq(0, 2, .1), col='lightgrey')
  plotQQ(getPlottingValues(ppCongruent, dataCongruent))
  plotQQ(getPlottingValues(ppIncongruent, dataIncongruent), colModel=3, colData=3)
  plotQQ(getPlottingValues(ppNeutral, dataNeutral), colModel=2, colData=2)
  legend('topright', c('C', 'N', 'I'), col=c(1,2,3), pch=c(15,15,15), lwd=2, bty='n', x.intersp = 0, lty=NA)
}

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}
