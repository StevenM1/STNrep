# Model 3 -----------------------------------------------------------------
load_model ("DDM","ddm.R") 
## drift and start point bias
transform.dmc <- function(par.df) 
{
  par.df$z = par.df$z + par.df$z_bias*par.df$z_sign
  par.df$v = par.df$v + par.df$v_bias*par.df$z_sign
  
  par.df$d <- par.df$d*par.df$t0  # proportional to t0, bounded -1 to 1
  par.df[,c("a","v","t0","z","d","sz","sv","st0")]
}

# Find bias trials
is_bias <- empty.map(list(S=c('stimleft', 'stimright'), coherence=c("low", "high"),
                          cue=c('cueleft', 'cueright', 'cueneutral'), R=c('r1', 'r2')),
                     levels=c('biastrial', 'neutraltrial'))
is_bias[grepl('cueleft', names(is_bias))] = 'biastrial'
is_bias[grepl('cueright', names(is_bias))] = 'biastrial'
is_bias <- assign.map(map=is_bias, value='neutraltrial')

# Find valid/congruent trials
is_valid <- empty.map(list(S=c('stimleft', 'stimright'), coherence=c("low", "high"),
                           cue=c('cueleft', 'cueright', 'cueneutral'), R=c('r1', 'r2')),
                      levels=c('Cvalid', 'Cinvalid', 'Cneutral'))
is_valid[grepl('cueleft', names(is_valid)) & grepl('stimleft', names(is_valid))] = 'Cvalid'
is_valid[grepl('cueleft', names(is_valid)) & grepl('stimright', names(is_valid))] = 'Cinvalid'
is_valid[grepl('cueright', names(is_valid)) & grepl('stimleft', names(is_valid))] = 'Cinvalid'
is_valid[grepl('cueright', names(is_valid)) & grepl('stimright', names(is_valid))] = 'Cvalid'
is_valid <- assign.map(map=is_valid, value='Cneutral')


p.map <- list(a="1", v="coherence", v_bias = "IS_BIAS", z="1",
              z_bias="IS_BIAS", z_sign="IS_VALID", d="1",sz="1",sv="1",t0="1",st0="1")
model <- model.dmc(
  p.map     = p.map,
  match.map = list(M=list(stimleft="r1", stimright="r2"),
                   IS_BIAS=is_bias, IS_VALID=is_valid),
  factors=list(S=c("stimleft", "stimright"),
               coherence=c("low", "high"),
               cue=c('cueleft', 'cueright', 'cueneutral')),
  constants = c(st0=0,d=0, sz=0, sv=0,
                z_sign.Cvalid=1, z_sign.Cinvalid=-1, z_sign.Cneutral=0,
                z_bias.neutraltrial=0, v_bias.neutraltrial=0),
  responses = c("r1","r2"),
  type = "rd")


# Prior distributions
p.mean <- c(a=2, v.low=3, v.high=3, v_bias.biastrial=0, z=0.5, z_bias.biastrial=0.0, t0=0.3)
p.scale <-c(a=2, v.low=3, v.high=3, v_bias.biastrial=2, z=0.1, z_bias.biastrial=0.1, t0=0.3)
p.prior <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1=p.mean,
  p2=p.scale,
  lower=c(0,-5, -5, -5, .3, -.3, 0),
  upper=c(5, 7,  7,  7, .7,  .3, 1)
)
