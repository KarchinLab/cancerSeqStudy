# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'mcores', 'c', 1, 'integer',
    'output', 'o', 1, 'character',
    'ratioMetric', 'r', 2, 'double',
    'driverFrac', 'd', 2, 'double',
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

#############################
# Analyze power and false positives
# when using a beta-binomial model
#############################
smgBbdFullAnalysis <- function(mu, cv, Leff, signif.level, effect.size, 
                               desired.power, samp.sizes){

  # find the power and numer of samples needed for a desired power
  powerResult <- smgBbdRequiredSampleSize(desired.power, mu, cv, samp.sizes, 
                                          effect.size, signif.level, Leff)
  bbd.samp.size.min <- powerResult$samp.size.min
  bbd.samp.size.max <- powerResult$samp.size.max
  power.result.bbd <- powerResult$power
  
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(mu, cv)
  
  # find expected number of false positives
  fp.result <- smg.binom.false.pos(params$alpha, params$beta, samp.sizes, Leff, 
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


ratiometricBbdFullAnalysis <- function(p, cv, mu, Leff, signif.level, effect.size, 
                                       desired.power, samp.sizes, Df){

  # find the power and numer of samples needed for a desired power
  powerResult <- ratiometricBbdRequiredSampleSize(p, cv, desired.power, samp.sizes, mu,
                                                  effect.size, Df, signif.level, Leff)
  bbd.samp.size.min <- powerResult$samp.size.min
  bbd.samp.size.max <- powerResult$samp.size.max
  power.result.bbd <- powerResult$power
  
  # get alpha and beta parameterization
  # for beta-binomial
  params <- rateCvToAlphaBeta(p, cv)
  
  # find expected number of false positives
  fp.result <- ratiometric.binom.false.pos(params$alpha, params$beta, samp.sizes, 
                                           mu, Leff, signif.level=signif.level)
  
  # save binomial data
  tmp.df <- data.frame(sample.size=samp.sizes)
  tmp.df["Power"] <- power.result.bbd
  tmp.df['sample min'] <- bbd.samp.size.min
  tmp.df['sample max'] <- bbd.samp.size.max
  tmp.df['CV'] <- cv
  tmp.df['signif.level'] <- signif.level
  tmp.df['effect.size'] <- effect.size
  tmp.df['mutation.rate'] <- mu
  tmp.df["Ratio-metric P"] <- p
  tmp.df["Driver fraction"] <- Df
  tmp.df["FP"] <- fp.result
  
  return(tmp.df)
}

######################
# Estimate power with a binomial model
######################

smgBinomFullAnalysis <- function(mu, Leff, signif.level, effect.size, 
                                 desired.power, samp.sizes){
  # calculate power
  power.result.binom <- smg.binom.power(mu, samp.sizes, Leff, 
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

ratiometricBinomFullAnalysis <- function(p, mu, Leff, signif.level, effect.size, 
                                         desired.power, samp.sizes, Df){
  # calculate power
  power.result.binom <- ratiometric.binom.power(p, samp.sizes, mu, Leff=Leff, 
                                                Df=Df, signif.level=signif.level,
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
  tmp.df["Ratio-metric P"] <- p
  tmp.df["Driver fraction"] <- Df
  tmp.df["FP"] <- NA
  
  return(tmp.df)
}

#############################
# run the analysis
#############################

runSmgAnalysisList <- function(x, samp.sizes, 
                               desired.power=.9, Leff=1500*3/4, 
                               possible.cvs=c()){
  # unpack the parameters
  mymu <- x[1]
  myeffect.size <- x[2]
  myalpha.level <- x[3]
  
  # run analysis
  result.df <- runSmgAnalysis(mymu, myeffect.size, myalpha.level,
                              samp.sizes, desired.power, Leff, possible.cvs)

  return(result.df)
}

runRatiometricAnalysisList <- function(x, samp.sizes, 
                                       Df=1.0, desired.power=.9, Leff=1500*3/4, 
                                       possible.cvs=c()){
  # unpack the parameters
  myp <- x[1]
  mymu <- x[2]
  myeffect.size <- x[3]
  myalpha.level <- x[4]
  
  # run analysis
  result.df <- runRatiometricAnalysis(myp, mymu, myeffect.size, myalpha.level,
                                      samp.sizes, desired.power, Leff, Df, possible.cvs)

  return(result.df)
}

#' This function unpacks a vector x which contains many combinations of the mutation
#' rate, effect.size, and significance level. The purpose of this function is parallelized
#' code running over a list of parameters. If you are not parallelizing, then use the 
#' runAnalysis function.
runAnalysisList <- function(x, analysisType="smg", Leff=1500*3/4, Df=1.0, ...){ 
  if (analysisType=="smg"){
    result <- runSmgAnalysisList(x, Leff=Leff, ...) 
  } else{
    result <- runRatiometricAnalysisList(x, Leff=Leff, Df=Df, ...) 
  }
  return(result)
}

#' Runs the entire power and false positive analysis pipeline.
runSmgAnalysis <- function(mu, effect.size, signif.level,
                           samp.sizes, desired.power=.9, 
                           Leff=1500*3/4, possible.cvs=c()){
  # run beta-binomial model
  result.df <- data.frame()
  for (mycv in possible.cvs){   
    # calculate false positives and power
    tmp.df <- smgBbdFullAnalysis(mu, mycv, Leff, signif.level, effect.size, 
                                 desired.power, samp.sizes)
    result.df <- rbind(result.df, tmp.df)
  }
  
  # save binomial data
  tmp.df <- smgBinomFullAnalysis(mu, Leff, signif.level, effect.size, desired.power, samp.sizes)
  result.df <- rbind(result.df, tmp.df)
  
  return(result.df)
}

#' Runs the entire power and false positive analysis pipeline.
runRatiometricAnalysis <- function(p, mu, effect.size, signif.level,
                                   samp.sizes, desired.power=.9, 
                                   Leff=1500*3/4, Df=1.0, possible.cvs=c()){
  # run beta-binomial model
  result.df <- data.frame()
  for (mycv in possible.cvs){   
    # calculate false positives and power
    tmp.df <- ratiometricBbdFullAnalysis(p, mycv, mu, Leff, signif.level, effect.size, 
                                         desired.power, samp.sizes, Df)
    result.df <- rbind(result.df, tmp.df)
  }
  
  # save binomial data
  tmp.df <- ratiometricBinomFullAnalysis(p, mu, Leff, signif.level, 
                                         effect.size, desired.power, 
                                         samp.sizes, Df)
  result.df <- rbind(result.df, tmp.df)
  
  return(result.df)
}

#' Utility function to get location of script when
#' executing Rscript.
thisfile <- function() {
    cmdArgs = commandArgs(trailingOnly = FALSE)
    needle = "--file="
    match = grep(needle, cmdArgs)
    if (length(match) > 0) {
        # Rscript
        return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
        ls_vars = ls(sys.frames()[[1]])
        if ("fileName" %in% ls_vars) {
            # Source'd via RStudio
            return(normalizePath(sys.frames()[[1]]$fileName)) 
        } else {
            # Source'd via R console
            return(normalizePath(sys.frames()[[1]]$ofile))
        }
    }
}

# Run as a script if arguments provided
if (!is.null(opt$ARGS)){
  # load other R files
  script.path <- thisfile()
  script.dir <- dirname(script.path)
  source(file.path(script.dir, "smg.R"))
  source(file.path(script.dir, "ratioMetric.R"))
  
  #############################
  # define the model params
  #############################
  # figure out whether to run a ratio-metric analysis
  # or mutation rate analysis
  if (is.null(opt$ratioMetric)){
    cmdType <- "smg" 
  } else {
    cmdType <- "ratio-metric" 
  }
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
  alpha.levels <- c(5e-6)  # list for level of significance

  # fraction of driver mutations being a certain
  # category of interest for ratio-metric methods
  driver.frac <- opt$driverFrac
  
  # setting up the sample sizes to check
  N <- 25000
  by.step <- 10
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
        if(cmdType=="smg"){
          param.list[[counter]] <- c(rate[i], effect.size, alpha.level)
        }else {
          param.list[[counter]] <- c(opt$ratioMetric, rate[i], effect.size, alpha.level)
        }
        counter <- counter + 1
      }
    }
  }
    
  ############################
  # run analysis
  ############################
  result.list <- mclapply(param.list, runAnalysisList, mc.cores=opt$mcores,
                          analysisType=cmdType, samp.sizes=samp.sizes, 
                          desired.power=desired.power, Df=driver.frac,
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
