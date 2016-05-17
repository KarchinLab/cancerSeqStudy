suppressPackageStartupMessages(library(VGAM))

##################################
# Power calculatations
##################################

#' calculates the power in a binomial power model
#' for ratio-metric approach
#'
#' @param p background proportion of total mutations falling into specific category
#' @param N vector of sample sizes
#' @param mu per base rate of mutation 
#' @param Df fraction of driver mutations that are the specific one of interest
#' @param Leff gene length in bases
#' @param r effect size for power analysis
#' @param signif.level alpha level for power analysis
#' @return vector containing power for each sample size
ratiometric.binom.power <- function(p, N, mu, 
                                    Df=1.0, Leff=1500*3/4, r=.02,
                                    signif.level=5e-6){
  # figure out the target mutation rate for effect size is
  muEffect <- 1 - ((1-mu)^(Leff) - r)^(1/Leff)
  # Calculate the discrepancy between the background and
  # target effect size
  muDiff <- muEffect - mu
  # given the mutation rates calculate the target effect
  # size for a ratio-metric method
  pEffect <- (mu*p + Df*muDiff) / muEffect
  
  # iterate over the number of samples
  power <- c()
  for(i in N){
    # step one, find the # of mutations where
    # it is expected to occur at least 90% of the time
    j <- 1
    while(j){
      prob <- pbinom(j-1, Leff*i, muEffect)
      if(prob >= .1){
        mutEff <- j
        break
      }
      j <- j+1
    }
    
    # step two, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbinom(j-1, mutEff, p)
      if(pval <= signif.level){
        Xc <- j
        break
      }
      j <- j+1
    }
    
    # step three, calculate power
    prob <- 1-pbinom(Xc-1, mutEff, pEffect)
    power <- c(power, prob) 
  }
  return(power)
}

#' Calculates the power in a ratio-metric approach using 
#' a beta-binomial power model.
#' 
#' The alpha and beta parameterize a proportion out of the
#' total mutations in a gene, rather than a mutation rate per base.
#'
#' @param my.alpha alpha parameter for beta binomial
#' @param my.beta beta parameter for beta binomial
#' @param N vector of sample sizes
#' @param mu per base rate of mutation 
#' @param Leff length of gene in bases
#' @param r effect size for power analysis
#' @param signif.level alpha level for power analysis
#' @return vector containing power for each sample size
ratiometric.bbd.power <- function(my.alpha, my.beta, 
                                  N, mu, Df=1.0,
                                  Leff=1500*3/4, r=.02,
                                  signif.level=5e-6){
  # figure out what the ratio-metric probability is from
  # the alpha and beta parameters
  p <- my.alpha / (my.alpha + my.beta)
  # figure out the target mutation rate for effect size is
  muEffect <- 1 - ((1-mu)^(Leff) - r)^(1/Leff)
  # Calculate the discrepancy between the background and
  # target effect size
  muDiff <- muEffect - mu
  # given the mutation rates calculate the target effect
  # size for a ratio-metric method
  pEffect <- (mu*p + Df*muDiff) / muEffect
  
  # iterate over the number of samples
  power <- c()
  for(i in N){
    # step one, find the # of mutations where
    # it is expected to occur at least 90% of the time
    j <- 1
    while(j){
      prob <- pbinom(j-1, Leff*i, muEffect)
      if(prob >= .1){
        mutEff <- j
        break
      }
      j <- j+1
    }
    
    # step two, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbetabinom.ab(j-1, mutEff, my.alpha, my.beta)
      if(pval <= signif.level){
        Xc <- j
        break
      }
      j <- j+1
    }
    
    # step three, calculate power
    prob <- 1-pbinom(Xc-1, mutEff, pEffect)
    power <- c(power, prob) 
  }
  return(power)
}

##################################
# Estimated false postives
##################################

#' calculates the false positives for a binomial model of
#' a ratio-metric feature.
#'
#' @param my.alpha alpha parameter for beta binomial
#' @param my.beta beta parameter for beta binomial
#' @param N vector of # samples to calculate power for
#' @param mu mutation rate per base
#' @param L gene CDS length in bases
#' @param num.genes number of genes that are tested
#' @param signif.level alpha level for power analysis
ratiometric.binom.false.pos <- function(my.alpha, my.beta,
                                        N, mu, Leff=1500*3/4,
                                        num.genes=18500,
                                        signif.level=5e-6){
  # calculate the ratio-metric fraction from alpha and beta
  p <- my.alpha / (my.alpha + my.beta)
  # examine power of binomial test
  # first find critical value based on binomial distribution
  power <- c()
  falsePositives <- c()
  for(i in N){
    # step one, find the # of mutations where
    # it is expected to occur at least 90% of the time
    #j <- 1
    #while(j){
    #  prob <- pbinom(j-1, L*i, mu)
    #  if(prob >= .1){
    #    mutEff <- j
    #    break
    #  }
    #  j <- j+1
    #}
    mutEff <- ceiling(Leff*i*mu)
    
    # step one, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbinom(j-1, mutEff, p)
      if(pval <= signif.level){
        Xc <- j
        break
      }
      j <- j+1
    }
    
    # step two, calculate false positives if overdispersion
    fp <- 1 - pbetabinom.ab(Xc-1, mutEff, my.alpha, my.beta)
    falsePositives <- c(falsePositives, num.genes*fp)
  }
  return(falsePositives)
}

############################
# Calculate required samples size
############################

#' Calculates the smallest sample size to detect driver genes for which
#' there is sufficient power using a binomial model for ratio-metric features.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a non-silent
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param p the background fraction of total mutations represented by the ratio-metric feature (e.g. inactivating mutations / total)
#' @param desired.power A floating point number indicating desired power
#' @param possible.samp.sizes vector of possible number of cancer samples in study
#' @param mu mutation rate per base
#' @param effect.size fraction of samples above background mutation rate
#' @param signif.level significance level for binomial test
#' @param L gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
ratiometricBinomRequiredSampleSize <- function(p, desired.power, possible.samp.sizes, mu,
                                               effect.size, Df=1.0, signif.lvl=5e-6, Leff=1500*3/4){
  # calculate power
  power.result.ratio <- ratiometric.binom.power(p, possible.samp.sizes, mu, Leff, 
                                                Df=Df, signif.level=signif.lvl,
                                                r=effect.size)
  ratiometric.samp.size.min <- possible.samp.sizes[min(which(power.result.ratio>=desired.power))]
  ratiometric.samp.size.max <- possible.samp.sizes[max(which(power.result.ratio<desired.power))+1]
  
  # return result
  result <- list(samp.size.min=ratiometric.samp.size.min, samp.size.max=ratiometric.samp.size.max,
                 power=power.result.ratio, sample.sizes=possible.samp.sizes)
  return(result)
}

#' Calculates the smallest sample size to detect driver genes for which
#' there is sufficient power using a beta-binomial model for ratio-metric features.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a non-silent
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param p the background fraction of total mutations represented by the ratio-metric feature (e.g. inactivating mutations / total)
#' @param cv the coefficient of variation for the parameter p
#' @param desired.power A floating point number indicating desired power
#' @param possible.samp.sizes vector of possible number of cancer samples in study
#' @param mu mutation rate per base
#' @param effect.size fraction of samples above background mutation rate
#' @param signif.level significance level for binomial test
#' @param L gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
ratiometricBbdRequiredSampleSize <- function(p, cv, desired.power, possible.samp.sizes, mu,
                                             effect.size, Df=1.0, signif.lvl=5e-6, Leff=1500*3/4){
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(p, cv)
  
  # calculate power
  power.result.ratio <- ratiometric.bbd.power(params$alpha, params$beta, 
                                              possible.samp.sizes,
                                              mu, Leff, Df=Df,
                                              signif.level=signif.lvl,
                                              r=effect.size)
  ratiometric.samp.size.min <- possible.samp.sizes[min(which(power.result.ratio>=desired.power))]
  ratiometric.samp.size.max <- possible.samp.sizes[max(which(power.result.ratio<desired.power))+1]
  
  # return result
  result <- list(samp.size.min=ratiometric.samp.size.min, samp.size.max=ratiometric.samp.size.max,
                 power=power.result.ratio, sample.sizes=possible.samp.sizes)
  return(result)
}

####################################
# Functions to calculate the minimum effect size
# with a given power
#####################################

#' Calculates the effect size of a driver gene according to a binomial model
#' of ratio-metric features for which there is sufficient power.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a 
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param possible.effect.sizes vector of effect sizes
#' @param desired.power A floating point number indicating desired power
#' @param p the background fraction of total mutations represented by the ratio-metric feature (e.g. inactivating mutations / total)
#' @param mu Mutation rate per base
#' @param samp.size number of cancer samples in study
#' @param signif.level significance level for binomial test
#' @param Leff effective gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
ratiometricBinomPoweredEffectSize <- function(possible.effect.sizes, desired.power, p, mu, 
                                              samp.size, Df=1.0, signif.level=5e-6, Leff=1500*3/4) {
  # calculate the power for each effect size
  pow.vec <- c()
  for(effect.size in possible.effect.sizes){
    pow <- ratiometric.binom.power(p, samp.size, mu, Leff, 
                                   Df=Df, signif.level=signif.level,
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

#' Calculates the effect size of a driver gene according to a beta-binomial model
#' of ratio-metric features for which there is sufficient power.
#' 
#' Effect size is measures as the fraction of sample/patient cancers with a 
#' mutation in a driver gene above the background mutation rate.
#' 
#' @param possible.effect.sizes vector of effect sizes
#' @param desired.power A floating point number indicating desired power
#' @param p the background fraction of total mutations represented by the ratio-metric feature (e.g. inactivating mutations / total)
#' @param cv the coefficient of variation for the parameter p
#' @param mu Mutation rate per base
#' @param samp.size number of cancer samples in study
#' @param signif.level significance level for binomial test
#' @param Leff effective gene length of CDS in bases for an average gene
#' @return List containing the smallest effect size with sufficient power
ratiometricBbdPoweredEffectSize <- function(possible.effect.sizes, desired.power, p, cv, mu, 
                                            samp.size, Df=1.0, signif.level=5e-6, Leff=1500*3/4) {
  # figure out alpha/beta for beta-binomial
  params <- rateCvToAlphaBeta(p, cv)
  
  # calculate the power for each effect size
  pow.vec <- c()
  for(effect.size in possible.effect.sizes){
    pow <- ratiometric.bbd.power(params$alpha, params$beta, 
                                 samp.size, mu, Leff, 
                                 Df=Df, signif.level=signif.level,
                                 r=effect.size)
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

