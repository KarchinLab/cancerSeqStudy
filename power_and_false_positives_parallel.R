library(VGAM)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(parallel)

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'mcores', 'c', 1, 'integer',
    'output', 'o', 1, 'character'
  ), byrow=TRUE, ncol=4)
  opt = getopt(spec)
  # print out help msg
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  }
} else {
  opt <- list(ARGS=NULL)
}

################
# calculates the power in a binomial power model
#
# Parameters
# my.mu : per base rate of mutation for binomial
# N : maximum number of sample to calculate power for
# r : effect size for power analysis
# alphaLevel : alpha level for power analysis
# by : discretization of sample sizes
################
binom.power <- function(my.mu,
                        Leff,
                        N=2000,
                        r=.02,
                        alpha.level=1e-4,
                        by=10){
  # examine power of binomial test
  # first find critical value based on binomial distribution
  # Calculate power for various sizes with different effects
  muEffect <- 1 - ((1-my.mu)^Leff - r)^(1/Leff)
  power <- c()
  falsePositives <- c()
  for(i in seq(by, N, by=by)){
    # step one, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbinom(j, Leff*i, my.mu)
      if(pval < alpha.level){
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

################
# calculates the false positives in a binomial model
# if there is over-diserspion
#
# Parameters
# my.mu : per base rate of mutation for binomial
# my.alpha : alpha parameter for beta binomial
# my.beta : beta parameter for beta binomial
# Leff : effective gene length in bases
# num.genes : number of genes that are tested
# N : maximum number of sample to calculate power for
# r : effect size for power analysis
# alphaLevel : alpha level for power analysis
# by : discretization of sample sizes
################
binom.false.pos <- function(my.mu, my.alpha, my.beta,
                            Leff,
                            num.genes=18500,
                            N=2000,
                            r=.02,
                            alpha.level=1e-4,
                            by=10){
  # examine power of binomial test
  # first find critical value based on binomial distribution
  # Calculate power for various sizes with different effects
  #alphaLevel <- 1e-4  # alpha level for genome-wide significance
  #N <- 2000  # potential maximum number of samples for study
  #r <- .02  # mutated gene per patient elevation above background
  muEffect <- 1 - ((1-my.mu)^Leff - r)^(1/Leff)
  power <- c()
  falsePositives <- c()
  for(i in seq(by, N, by=by)){
    # step one, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbinom(j, Leff*i, my.mu)
      if(pval < alpha.level){
        Xc <- j
        break
      }
      j <- j+1
    }
    
    # step two, calculate power
    p <- 1-pbinom(Xc-1, Leff*i, muEffect)
    power <- c(power, p)
    
    # step three, calculate false positives if overdispersion
    fp <- 1 - pbetabinom.ab(Xc-1, Leff*i, my.alpha, my.beta)
    falsePositives <- c(falsePositives, num.genes*fp)
  }
  return(falsePositives)
}

################
# calculates the power in a beta-binomial model
#
# Parameters
# my.mu : per base rate of mutation for binomial
# my.alpha : alpha parameter for beta binomial
# my.beta : beta parameter for beta binomial
# Leff : effective gene length in bases
# num.genes : number of genes that are tested
# N : maximum number of sample to calculate power for
# r : effect size for power analysis
# alphaLevel : alpha level for power analysis
# by : discretization of sample sizes
################
bbd.power <- function(my.mu, my.alpha, my.beta,
                      Leff,
                      num.genes=18500,
                      N=2000,
                      r=.02,
                      alpha.level=1e-4,
                      by=10){
  # examine power of binomial test
  # first find critical value based on binomial distribution
  # Calculate power for various sizes with different effects
  #alphaLevel <- 1e-4  # alpha level for genome-wide significance
  #N <- 2000  # potential maximum number of samples for study
  #r <- .02  # mutated gene per patient elevation above background
  muEffect <- 1 - ((1-my.mu)^Leff - r)^(1/Leff)
  print(muEffect)
  print(my.mu)
  power <- c()
  falsePositives <- c()
  for(i in seq(by, N, by=by)){
    # step one, find critical threshold
    j <- 1
    while(j){
      pval <- 1-pbetabinom.ab(j, Leff*i, my.alpha, my.beta)
      if(pval < alpha.level){
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
# run the analysis
#############################
run.analysis <- function(x){
  result.df <- data.frame()
  
  for (mycv in possible.cv){
    # unpack the parameters
    mypi <- x[1]
    myeffect.size <- x[2]
    myalpha.level <- x[3]
    
    # define beta binomial params
    ab <- mypi * (1-mypi) / (mycv*mypi)^2 - 1
    my.alpha <- mypi * ab
    my.beta <- (1-mypi)*ab
    my.mu <- mypi
    
    # choose which mutation rate to use
    #my.mu <- mypi[i]
    #my.alpha <- alpha[i]
    #my.beta <- beta[i]
    
    # calc power
    #power.result.binom <- binom.power(my.mu, Leff, alpha.level=alpha.level, N=N, r=effect.size)
    power.result.bbd <- bbd.power(my.mu, my.alpha, my.beta, Leff, alpha.level=myalpha.level, 
                                  by=by.step, N=N, r=myeffect.size)
    
    # find min/max samples to achieve desired power
    #binom.samp.size.min <- samp.sizes[min(which(power.result.binom>=desired.power))]
    #binom.samp.size.max <- samp.sizes[max(which(power.result.binom<desired.power))+1]
    bbd.samp.size.min <- samp.sizes[min(which(power.result.bbd>=desired.power))]
    bbd.samp.size.max <- samp.sizes[max(which(power.result.bbd<desired.power))+1]
    
    # find expected number of false positives
    fp.result <- binom.false.pos(my.mu, my.alpha, my.beta, Leff, alpha.level=myalpha.level, 
                                 by=by.step, N=N, r=myeffect.size)
    
    # save binomial data
    tmp.df <- data.frame(sample.size=samp.sizes)
    tmp.df["Power"] <- power.result.bbd
    tmp.df['sample min'] <- bbd.samp.size.min
    tmp.df['sample max'] <- bbd.samp.size.max
    tmp.df['CV'] <- mycv
    tmp.df['alpha.level'] <- myalpha.level
    tmp.df['effect.size'] <- myeffect.size
    tmp.df['mutation.rate'] <- my.mu
    tmp.df["FP"] <- fp.result
    result.df <- rbind(result.df, tmp.df)
  }
  
  # save binomial data
  power.result.binom <- binom.power(my.mu, Leff, alpha.level=myalpha.level,
                                    by=by.step, N=N, r=myeffect.size)
  binom.samp.size.min <- samp.sizes[min(which(power.result.binom>=desired.power))]
  binom.samp.size.max <- samp.sizes[max(which(power.result.binom<desired.power))+1]
  tmp.df["Power"] <- power.result.binom
  tmp.df['sample min'] <- binom.samp.size.min
  tmp.df['sample max'] <- binom.samp.size.max
  tmp.df['CV'] <- 0
  result.df <- rbind(result.df, tmp.df)
  return(result.df)
}

#############################
# define the model params
#############################
rate <<- c(.1e-6, .2e-6, .3e-6, .5e-6, 7e-6, 1e-6, 1.5e-6, 2e-6, 3e-6, 
           4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 10e-6)
fg <<- 3.9  # an adjustment factor that lawrence et al used for variable gene length
rate <<- fg*rate

nonsilentFactor <<- 3/4
L <<- 1500  # same length as used in paper
Leff <<- L * nonsilentFactor
N <<- 20000
by.step <<- 25
samp.sizes <<- seq(by.step, N, by=by.step)
desired.power <<- .9
possible.cv <<- c(.05, .1, .2)
effect.sizes <<- c(.01, .02, .05)
alpha.levels <<- c(1e-4, 5e-6)

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
      param.list[[counter]] <- c(rate[i], effect.size, alpha.level, possible.cv)
      counter <- counter + 1
    }
  }
}

############################
# run analysis
############################
result.list <- mclapply(param.list, run.analysis, mc.cores=opt$mcores)
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