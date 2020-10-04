# Model 2 -----------------------------------------------------------------
load_model ("DDM","ddm.R") 
## cue effect on drift only
transform.dmc <- function(par.df) 
{
  par.df$z = par.df$z + par.df$z_bias*par.df$z_sign
  par.df$v = par.df$v + par.df$v_bias*par.df$z_sign
  
  par.df$d <- par.df$d*par.df$t0  # proportional to t0, bounded -1 to 1
  par.df[,c("a","v","t0","z","d","sz","sv","st0")]
}

p.map <- list(a="1", v="coherence", v_bias = "IS_BIAS", z="1", 
              z_bias="IS_BIAS", z_sign="IS_VALID", d="1",sz="1",sv="1",t0="1",st0="1")
model <- model.dmc(
  p.map     = p.map,
  match.map = list(M=list(stimleft="r1", stimright="r2"),
                   IS_BIAS=is_bias, IS_VALID=is_valid),
  factors=list(S=c("stimleft", "stimright"), 
               coherence=c("low", "high"),               
               cue=c('cueleft', 'cueright', 'cueneutral')),
  constants = c(st0=0,d=0, #sz=0, sv=0,
                z_bias.biastrial=0,
                z_sign.Cvalid=1, z_sign.Cinvalid=-1, z_sign.Cneutral=0, 
                z_bias.neutraltrial=0, v_bias.neutraltrial=0),
  responses = c("r1","r2"),
  type = "rd")

# Prior distributions
p.mean <- c(a=2, v.low=3, v.high=3, v_bias.biastrial=0, z=0.5, t0=0.3, sv=0.1, sz=0.1)
p.scale <-c(a=2, v.low=3, v.high=3, v_bias.biastrial=1, z=0.1, t0=0.3, sv=0.5, sz=0.1)
p.prior <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p.mean,                           
  p2=p.scale,
  lower=c(0,-5, -5, -5, .3, 0, 0, 0),
  upper=c(5, 7,  7,  7, .7, 1, 5, 0.3)
)

# simulate model 2 to check if cue is correctly affecting v
# p.vector <- c("a"=1, "v.low"=1, "v.high"=1, "z"=0.5, "z_bias.biastrial"=0, "v_bias.biastrial"=1, "t0"=0, sv=0.1, sz=0.1)
# raw.data <- simulate.dmc(model, p.vector=p.vector, n = 500)
# par(mfrow=c(2,3))
# plot.cell.density(data.cell=raw.data[raw.data$S=="stimleft" & raw.data$cue=='cueneutral',],C="r1")
# plot.cell.density(data.cell=raw.data[raw.data$S=="stimleft" & raw.data$cue=='cueleft',],C="r1")
# plot.cell.density(data.cell=raw.data[raw.data$S=="stimleft" & raw.data$cue=='cueright',],C="r1")
# 
# plot.cell.density(data.cell=raw.data[raw.data$S=="stimright" & raw.data$cue=='cueneutral',],C="r2")
# plot.cell.density(data.cell=raw.data[raw.data$S=="stimright" & raw.data$cue=='cueleft',],C="r2")
# plot.cell.density(data.cell=raw.data[raw.data$S=="stimright" & raw.data$cue=='cueright',],C="r2")
