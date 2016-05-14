suppressPackageStartupMessages(library(VGAM))

##################################
# Power calculatations
##################################

#' calculates the power in a binomial power model
#' for significantly mutated genes
#'
#' @param my.mu per base rate of mutation for binomial
#' @param N vector of sample sizes
#' @param r effect size for power analysis
#' @param signif.level alpha level for power analysis
#' @return vector containing power for each sample size
smg.binom.power <- function(my.mu,
                            N,
                            Leff=1500*3/4,
                            r=.02,
                            signif.level=5e-6){
  # examine power of binomial test
  # first find critical value based on binomial distribution
  # Calculate power for various sizes with different effects
  muEffect <- 1 - ((1-my.mu)^Leff - r)^(1/Leff)
  power <- c()
  falsePositives <- c()
  for(i in N){
    # step one, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbinom(j-1, Leff*i, my.mu)
      if(pval <= signif.level){
        Xc <- j
        break
      }
      j <- j+1
    }
    
    # step two, calculate power
    p <- 1-pbinom(Xc-1, Leff*i, muEffect)
    power <- c(power, p)
    
  }

  return(power)
}

#' calculates the power in a beta-binomial model for
#' significantly mutated genes.
#'
#' @param my.alpha alpha parameter for beta binomial
#' @param my.beta beta parameter for beta binomial
#' @param N maximum number of sample to calculate power for
#' @param Leff effective gene length in bases
#' @param r effect size for power analysis
#' @param signif.level alpha level for power analysis
smg.bbd.power <- function(my.alpha, my.beta,
                          N,
                          Leff=1500*3/4,
                          r=.02,
                          signif.level=5e-6){
  # calc the mutation rate from alpha/beta
  my.mu <- my.alpha / (my.alpha + my.beta)
  # examine power of binomial test
  # first find critical value based on binomial distribution
  # Calculate power for various sizes with different effects
  muEffect <- 1 - ((1-my.mu)^Leff - r)^(1/Leff)
  power <- c()
  falsePositives <- c()
  for(i in N){
    # step one, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbetabinom.ab(j-1, Leff*i, my.alpha, my.beta)
      if(pval <= signif.level){
        Xc <- j
        break
      }
      j <- j+1
    }
    
    # step two, calculate power
    p <- 1-pbinom(Xc-1, Leff*i, muEffect)
    power <- c(power, p)

  }
  return(power)
}

###################################
# Estimating false positives
###################################

#' calculates the false positives in a binomial model
#' for identifying significantly mutated genes if 
#' there is over-diserspion.
#'
#' @param my.alpha alpha parameter for beta binomial
#' @param my.beta beta parameter for beta binomial
#' @param N vector of # samples to calculate power for
#' @param Leff effective gene length in bases
#' @param num.genes number of genes that are tested
#' @param signif.level alpha level for power analysis
smg.binom.false.pos <- function(my.alpha, my.beta,
                                N, Leff=1500*3/4,
                                num.genes=18500,
                                signif.level=5e-6){
  # calculate mutation rate from alpha/beta
  my.mu <- my.alpha / (my.alpha + my.beta)
  # examine power of binomial test
  # first find critical value based on binomial distribution
  power <- c()
  falsePositives <- c()
  for(i in N){
    # step one, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbinom(j-1, Leff*i, my.mu)
      if(pval <= signif.level){
        Xc <- j
        break
      }
      j <- j+1
    }
    
    # step two, calculate false positives if overdispersion
    fp <- 1 - pbetabinom.ab(Xc-1, Leff*i, my.alpha, my.beta)
    falsePositives <- c(falsePositives, num.genes*fp)
  }
  return(falsePositives)
}

############################
# Calculate required samples size
############################

#' Calculates the smallest sample size to detect driver genes for which
#' there is sufficient power using a beta-binomial model.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a non-silent
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param desired.power A floating point number indicating desired power
#' @param mu Mutation rate per base
#' @param cv Coefficient of Variation surrounding the uncertaintly in mutation rate
#' @param possible.samp.sizes vector of possible number of cancer samples in study
#' @param effect.size fraction of samples above background mutation rate
#' @param signif.level significance level for binomial test
#' @param Leff effective gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
smgBbdRequiredSampleSize <- function(desired.power, mu, cv, possible.samp.sizes, 
                                     effect.size, signif.level=5e-6, Leff=1500*3/4){
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(mu, cv)
  
  # calc power
  power.result.bbd <- smg.bbd.power(params$alpha, params$beta, possible.samp.sizes, Leff, 
                                    signif.level=signif.level, r=effect.size)
  
  # find min/max samples to achieve desired power
  bbd.samp.size.min <- possible.samp.sizes[min(which(power.result.bbd>=desired.power))]
  bbd.samp.size.max <- possible.samp.sizes[max(which(power.result.bbd<desired.power))+1]
  
  # return result
  result <- list(samp.size.min=bbd.samp.size.min, samp.size.max=bbd.samp.size.max,
                 power=power.result.bbd, sample.sizes=possible.samp.sizes)
  return(result)
}

#' Calculates the smallest sample size to detect driver genes for which
#' there is sufficient power using a binomial model for mutation rate.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a non-silent
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param desired.power A floating point number indicating desired power
#' @param mu Mutation rate per base
#' @param possible.samp.sizes vector of possible number of cancer samples in study
#' @param effect.size fraction of samples above background mutation rate
#' @param signif.level significance level for binomial test
#' @param Leff effective gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
smgBinomRequiredSampleSize <- function(desired.power, mu, possible.samp.sizes,
                                       effect.size, signif.level=5e-6, Leff=1500*3/4){
  # calculate power
  power.result.binom <- smg.binom.power(mu, possible.samp.sizes, Leff, 
                                        signif.level=signif.level,
                                        r=effect.size)
  binom.samp.size.min <- possible.samp.sizes[min(which(power.result.binom>=desired.power))]
  binom.samp.size.max <- possible.samp.sizes[max(which(power.result.binom<desired.power))+1]
  
  # return result
  result <- list(samp.size.min=binom.samp.size.min, samp.size.max=binom.samp.size.max,
                 power=power.result.binom, sample.sizes=possible.samp.sizes)
  return(result)
}


################################
# Calculates the effect size which has power
################################

#' Calculates the smallest effect size in a driver gene for which
#' there is sufficient power using a significantly mutated gene
#' approach with a beta-binomial model.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a non-silent
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param possible.effect.sizes vector of effect sizes
#' @param desired.power A floating point number indicating desired power
#' @param mu Mutation rate per base
#' @param cv Coefficient of Variation surrounding the uncertaintly in mutation rate
#' @param samp.size number of cancer samples in study
#' @param signif.level significance level for binomial test
#' @param Leff effective gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
smgBbdPoweredEffectSize <- function(possible.effect.sizes, desired.power, mu, cv, samp.size, 
                                    signif.level=5e-6, Leff=1500*3/4) {
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(mu, cv)
  
  # calculate the power for each effect size
  pow.vec <- c()
  for(effect.size in possible.effect.sizes){
    # calc power
    pow <- smg.bbd.power(params$alpha, params$beta, samp.size, Leff, 
                         signif.level=signif.level, r=effect.size)
    pow.vec <- c(pow.vec, pow)
  }
  
  # find the effect size
  bbd.eff.size.min <- possible.effect.sizes[min(which(pow.vec>=desired.power))]
  bbd.eff.size.max <- possible.effect.sizes[max(which(pow.vec<desired.power))+1]
  
  # return result
  result <- list(eff.size.min=bbd.eff.size.min, eff.size.max=bbd.eff.size.max,
                 power=pow.vec, eff.size=possible.effect.sizes)
  return(result)
}

#' Calculates the minimum effect size (with sufficient power) of a driver gene according to a binomial model 
#' for significantly mutated genes.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a non-silent
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param possible.effect.sizes vector of effect sizes
#' @param desired.power A floating point number indicating desired power
#' @param mu Mutation rate per base
#' @param samp.size number of cancer samples in study
#' @param signif.level significance level for binomial test
#' @param Leff effective gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
smgBinomPoweredEffectSize <- function(possible.effect.sizes, desired.power, mu, samp.size, 
                                      signif.level=5e-6, Leff=1500*3/4) {
  # calculate the power for each effect size
  pow.vec <- c()
  for(effect.size in possible.effect.sizes){
    pow <- smg.binom.power(mu, samp.size, Leff, 
                           signif.level=signif.level,
                           r=effect.size)
    pow.vec <- c(pow.vec, pow)
  }
  
  # find the effect size
  binom.eff.size.min <- possible.effect.sizes[min(which(pow.vec>=desired.power))]
  binom.eff.size.max <- possible.effect.sizes[max(which(pow.vec<desired.power))+1]
  
  # return result
  result <- list(eff.size.min=binom.eff.size.min, eff.size.max=binom.eff.size.max,
                 power=pow.vec, eff.size=possible.effect.sizes)
  return(result)
}
