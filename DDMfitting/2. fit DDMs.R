## Fit DDM
rm(list=ls())
setwd('~/bias/DDMfitting')
source ("dmc/dmc.R")
load_model ("DDM","ddm.R") 
source('utils.R')

# ## with sv, sz
# for(modelN in 1:3) {
#   # reload data
#   load('./behavior.RData')
#   
#   # load model
#   source(paste0('./model', modelN, 'sv.R'))
#   data.model <- data.model.dmc(dat, model)
#   
#   # sample
#   samples <- h.run.unstuck.dmc(h.samples.dmc(nmc = 500, p.prior, data.model, thin=1), p.migrate=0.05, cores=17, end.no.migrate = TRUE)
#   save(samples, file=paste0('samplesModel', modelN, 'sv.RData'))
#   samples <- h.samples.dmc(nmc=5000, add=FALSE, samples=samples)
#   samples <- h.RUN.dmc(samples, cores=17, cut.converge = 1.02)
#   save(samples, file=paste0('samplesModel', modelN, 'sv.RData'))
# }
# 
# ## No sv, sz
# for(modelN in 1:3) {
#   # reload data
#   load('./behavior.RData')
#   
#   # load model
#   source(paste0('./model', modelN, '.R'))
#   data.model <- data.model.dmc(dat, model)
#   
#   # sample
#   samples <- h.run.unstuck.dmc(h.samples.dmc(nmc = 500, p.prior, data.model, thin=1), p.migrate=0.05, cores=17, end.no.migrate = TRUE)
#   save(samples, file=paste0('samplesModel', modelN, '.RData'))
#   samples <- h.samples.dmc(nmc=5000, add=FALSE, samples=samples)
#   samples <- h.RUN.dmc(samples, cores=17, cut.converge = 1.02)
#   save(samples, file=paste0('samplesModel', modelN, '.RData'))
# }
# 
# 
# for(modelType in c('', 'sv')) {
#   for(modelN in 1:3) {
#     # reload data
#     print(modelN)
#     load('./behavior.RData')
#     
#     # load model
#     source(paste0('./model', modelN, modelType, '.R'))
#     load(paste0('samplesModel', modelN, modelType, '.RData'))
#     runDiagnostics(fn=paste0('model', modelN, modelType), dat=dat, samples=samples)
#   }
# }
# 


## with sv, sz, st0
for(modelN in 2:3) {
  # reload data
  load('./behavior.RData')
  
  # load model
  source(paste0('./model', modelN, 'st0.R'))
  data.model <- data.model.dmc(dat, model)
  
  # sample
  samples <- h.run.unstuck.dmc(h.samples.dmc(nmc = 500, p.prior, data.model, thin=1), p.migrate=0.05, cores=17, end.no.migrate = TRUE)
  save(samples, file=paste0('samplesModel', modelN, 'st0.RData'))
  samples <- h.samples.dmc(nmc=5000, add=FALSE, samples=samples)
  samples <- h.RUN.dmc(samples, cores=17, cut.converge = 1.02)
  save(samples, file=paste0('samplesModel', modelN, 'st0.RData'))
}


for(modelType in c('st0')) {
  for(modelN in 1:3) {
    # reload data
    print(modelN)
    load('./behavior.RData')
    
    # load model
    source(paste0('./model', modelN, modelType, '.R'))
    load(paste0('samplesModel', modelN, modelType, '.RData'))
    runDiagnostics(fn=paste0('model', modelN, modelType), dat=dat, samples=samples)
  }
}
