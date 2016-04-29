# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'mcores', 'c', 1, 'integer',
    'output', 'o', 1, 'character',
    'help', 'h', 0, 'logical'
  ), byrow=TRUE, ncol=4)
  opt = getopt(spec)
  # print out help msg
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  } else if (is.null(opt$mcores) | is.null(opt$output)){
    opt <- list(ARGS=NULL)
  }
} else {
  opt <- list(ARGS=NULL)
}

suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(parallel))

#' calculates the power in a binomial power model
#'
#' @param my.mu per base rate of mutation for binomial
#' @param N vector of sample sizes
#' @param r effect size for power analysis
#' @param signif.level alpha level for power analysis
#' @return vector containing power for each sample size
binom.power <- function(my.mu,
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
  #for(i in seq(by, N, by=by)){
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


#' calculates the false positives in a binomial model
#' if there is over-diserspion
#'
#' @param my.alpha alpha parameter for beta binomial
#' @param my.beta beta parameter for beta binomial
#' @param N vector of # samples to calculate power for
#' @param Leff effective gene length in bases
#' @param num.genes number of genes that are tested
#' @param signif.level alpha level for power analysis
binom.false.pos <- function(my.alpha, my.beta,
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


#' calculates the power in a beta-binomial model
#'
#' @param my.alpha alpha parameter for beta binomial
#' @param my.beta beta parameter for beta binomial
#' @param N maximum number of sample to calculate power for
#' @param Leff effective gene length in bases
#' @param r effect size for power analysis
#' @param signif.level alpha level for power analysis
bbd.power <- function(my.alpha, my.beta,
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

#############################
# Convert a rate and coefficient
# of variation parameter into
# the alpha and beta parameters
#############################

#' Converts mutation rate and coefficient of variation (CV) parameters
#' to equivalent alpha and beta parameters typically used for beta-binomial.
#' 
#' @param rate mutation rate
#' @param cv coefficient of variation for mutation rate
#' @return Param list containing alpha and beta
rateCvToAlphaBeta <- function(rate, cv) {
  ab <- rate * (1-rate) / (cv*rate)^2 - 1
  my.alpha <- rate * ab
  my.beta <- (1-rate)*ab
  return(list(alpha=my.alpha, beta=my.beta))
}

###################
# Functions to calculate the required sample size
##################

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
bbdRequiredSampleSize <- function(desired.power, mu, cv, possible.samp.sizes, 
                                  effect.size, signif.level=5e-6, Leff=1500*3/4){
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(mu, cv)
  
  # calc power
  power.result.bbd <- bbd.power(params$alpha, params$beta, possible.samp.sizes, Leff, 
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
#' there is sufficient power using a binomial model.
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
binomRequiredSampleSize <- function(desired.power, mu, possible.samp.sizes,
                                    effect.size, signif.level=5e-6, Leff=1500*3/4){
  # calculate power
  power.result.binom <- binom.power(mu, possible.samp.sizes, Leff, 
                                    signif.level=signif.level,
                                    r=effect.size)
  binom.samp.size.min <- possible.samp.sizes[min(which(power.result.binom>=desired.power))]
  binom.samp.size.max <- possible.samp.sizes[max(which(power.result.binom<desired.power))+1]
  
  # return result
  result <- list(samp.size.min=binom.samp.size.min, samp.size.max=binom.samp.size.max,
                 power=power.result.binom, sample.sizes=possible.samp.sizes)
  return(result)
}

####################################
# Functions to calculate the minimum effect size
# with a given power
#####################################

#' Calculates the smallest effect size in a driver gene for which
#' there is sufficient power using a beta-binomial model.
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
bbdPoweredEffectSize <- function(possible.effect.sizes, desired.power, mu, cv, samp.size, 
                                 signif.level=5e-6, Leff=1500*3/4) {
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(mu, cv)
  
  # calculate the power for each effect size
  pow.vec <- c()
  for(effect.size in possible.effect.sizes){
    # calc power
    pow <- bbd.power(params$alpha, params$beta, samp.size, Leff, 
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

#' Calculates the effect size of a driver gene according to a binomial for which
#' there is sufficient power.
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
binomPoweredEffectSize <- function(possible.effect.sizes, desired.power, mu, samp.size, 
                                   signif.level=5e-6, Leff=1500*3/4) {
  # calculate the power for each effect size
  pow.vec <- c()
  for(effect.size in possible.effect.sizes){
    pow <- binom.power(mu, samp.size, Leff, 
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

#############################
# Analyze power and false positives
# when using a beta-binomial model
#############################
bbdFullAnalysis <- function(mu, cv, Leff, signif.level, effect.size, 
                            desired.power, samp.sizes){

  # find the power and numer of samples needed for a desired power
  powerResult <- bbdRequiredSampleSize(desired.power, mu, cv, samp.sizes, 
                                       effect.size, signif.level, Leff)
  bbd.samp.size.min <- powerResult$samp.size.min
  bbd.samp.size.max <- powerResult$samp.size.max
  power.result.bbd <- powerResult$power
  
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(mu, cv)
  
  # find expected number of false positives
  fp.result <- binom.false.pos(params$alpha, params$beta, samp.sizes, Leff, 
                               signif.level=signif.level)
  
  # save binomial data
  tmp.df <- data.frame(sample.size=samp.sizes)
  tmp.df["Power"] <- power.result.bbd
  tmp.df['sample min'] <- bbd.samp.size.min
  tmp.df['sample max'] <- bbd.samp.size.max
  tmp.df['CV'] <- cv
  tmp.df['signif.level'] <- signif.level
  tmp.df['effect.size'] <- effect.size
  tmp.df['mutation.rate'] <- mu
  tmp.df["FP"] <- fp.result
  
  return(tmp.df)
}

binomFullAnalysis <- function(mu, Leff, signif.level, effect.size, 
                              desired.power, samp.sizes){
  # calculate power
  power.result.binom <- binom.power(mu, samp.sizes, Leff, 
                                    signif.level=signif.level,
                                    r=effect.size)
  binom.samp.size.min <- samp.sizes[min(which(power.result.binom>=desired.power))]
  binom.samp.size.max <- samp.sizes[max(which(power.result.binom<desired.power))+1]
  
  # record all power measurements
  tmp.df <- data.frame(sample.size=samp.sizes)
  tmp.df["Power"] <- power.result.binom
  tmp.df['sample min'] <- binom.samp.size.min
  tmp.df['sample max'] <- binom.samp.size.max
  tmp.df['CV'] <- 0
  tmp.df['signif.level'] <- signif.level
  tmp.df['effect.size'] <- effect.size
  tmp.df['mutation.rate'] <- mu
  tmp.df["FP"] <- NA
  
  return(tmp.df)
}

#############################
# run the analysis
#############################
#' This function unpacks a vector x which contains many combinations of the mutation
#' rate, effect.size, and significance level. The purpose of this function is parallelized
#' code running over a list of parameters. If you are not parallelizing, then use the 
#' runAnalysis function.
runAnalysisList <- function(x, samp.sizes, 
                            desired.power=.9, Leff=1500*3/4, 
                            possible.cvs=c()){
  # unpack the parameters
  mymu <- x[1]
  myeffect.size <- x[2]
  myalpha.level <- x[3]
  
  # run analysis
  result.df <- runAnalysis(mymu, myeffect.size, myalpha.level,
                           samp.sizes, desired.power, Leff, possible.cvs)

  return(result.df)
}

#' Runs the entire power and false positive analysis pipeline.
runAnalysis <- function(mu, effect.size, signif.level,
                        samp.sizes, desired.power=.9, 
                        Leff=1500*3/4, possible.cvs=c()){
  # run beta-binomial model
  result.df <- data.frame()
  for (mycv in possible.cvs){   
    # calculate false positives and power
    tmp.df <- bbdFullAnalysis(mu, mycv, Leff, signif.level, effect.size, 
                              desired.power, samp.sizes)
    result.df <- rbind(result.df, tmp.df)
  }
  
  # save binomial data
  tmp.df <- binomFullAnalysis(mu, Leff, signif.level, effect.size, desired.power, samp.sizes)
  result.df <- rbind(result.df, tmp.df)
  
  return(result.df)
}

# Run as a script if arguments provided
if (!is.null(opt$ARGS)){
  #############################
  # define the model params
  #############################
  # long list of rates to be evaluated
  rate <- c(.1e-6, .2e-6, .3e-6, .4e-6, .5e-6, .7e-6, .8e-6, 1e-6, 1.25e-6, 1.5e-6, 1.75e-6, 2e-6, 2.25e-6, 2.5e-6, 2.75e-6, 3e-6, 3.5e-6, 4e-6,
            4.5e-6, 5e-6, 5.5e-6, 6e-6, 6.5e-6, 7e-6, 7.5e-6, 8e-6, 8.5e-6, 9e-6, 10e-6, 11e-6, 12e-6)
  fg <- 3.9  # an adjustment factor that lawrence et al used for variable gene length
  rate <- fg*rate  # nominal rates are adjusted (will have to adjust back after analysis is done)
  
  # model parameters
  nonsilentFactor <- 3/4  # roughly the fraction 
  L <- 1500  # same length as used in lawrence et al. paper
  Leff <- L * nonsilentFactor
  desired.power <- .9  # aka 90% power
  possible.cvs <- c(.05, .1, .2)  # coefficient of variation for mutation rate per base
  effect.sizes <- c(.01, .02, .05)  # fraction of samples above background
  alpha.levels <- c(1e-4, 5e-6)  # level of significance
  
  # setting up the sample sizes to check
  N <- 25000
  by.step <- 25
  samp.sizes <- seq(by.step, N, by=by.step)  # grid of sample sizes to check
  
  ##################################
  # Loop through different params
  ##################################
  param.list <- list()
  counter <- 1
  for (i in 1:length(rate)){
    # loop over effect sizes
    for (effect.size in effect.sizes){
      # loop over alpha levels
      for (alpha.level in alpha.levels){
        param.list[[counter]] <- c(rate[i], effect.size, alpha.level)
        counter <- counter + 1
      }
    }
  }
  
  ############################
  # run analysis
  ############################
  result.list <- mclapply(param.list, runAnalysisList, mc.cores=opt$mcores,
                          samp.sizes=samp.sizes, desired.power=desired.power,
                          Leff=Leff, possible.cvs=possible.cvs)
  result.df <- do.call("rbind", result.list)
  
  # adjust mutation rates back to the average
  result.df$mutation.rate <- result.df$mutation.rate / fg
  # convert to factor
  result.df$mutation.rate <- factor(result.df$mutation.rate, levels=unique(result.df$mutation.rate))
  result.df$effect.size <- factor(result.df$effect.size, levels=unique(result.df$effect.size))
  
  ######################
  # Save result to text file
  ######################
  write.table(result.df, opt$output, sep='\t')
}
