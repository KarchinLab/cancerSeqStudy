library(VGAM)
library(reshape2)
library(ggplot2)
library(gridExtra)
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
# define the model params
#############################
pi <- c(.1e-6, .5e-6, 1e-6, 3e-6, 5e-6, 7e-6, 10e-6)
fg <- 3.9  # an adjustment factor that lawrence et al used for variable gene length
pi <- fg*pi

nonsilentFactor <- 3/4
L <- 1500  # same length as used in paper
Leff <- L * nonsilentFactor
N <- 10000
by=50
samp.sizes <- seq(by, N, by=by)
desired.power = .9
possible.cv = c(.05, .1, .2)
effect.sizes <- c(.01, .02, .05)
alpha.levels <- c(1e-4, 5e-6)


##################################
# Loop through different params
##################################
power.result <- data.frame()
false.pos.result <- data.frame()
result.df <- data.frame()
# loop over mutation rates
for (i in 1:length(pi)){
  # loop over effect sizes
  for (effect.size in effect.sizes){
    # loop over alpha levels
    for (alpha.level in alpha.levels){
      # loop over various coefficient of variations
      for (mycv in possible.cv){
        # define beta binomial params
        ab <- pi * (1-pi) / (mycv*pi)^2 - 1
        alpha <- pi * ab
        beta <- (1-pi)*ab
        
        # choose which mutation rate to use
        my.mu <- pi[i]
        my.alpha <- alpha[i]
        my.beta <- beta[i]
      
        # calc power
        #power.result.binom <- binom.power(my.mu, Leff, alpha.level=alpha.level, N=N, r=effect.size)
        power.result.bbd <- bbd.power(my.mu, my.alpha, my.beta, Leff, alpha.level=alpha.level, 
                                      by=by, N=N, r=effect.size)
        
        # find min/max samples to achieve desired power
        #binom.samp.size.min <- samp.sizes[min(which(power.result.binom>=desired.power))]
        #binom.samp.size.max <- samp.sizes[max(which(power.result.binom<desired.power))+1]
        bbd.samp.size.min <- samp.sizes[min(which(power.result.bbd>=desired.power))]
        bbd.samp.size.max <- samp.sizes[max(which(power.result.bbd<desired.power))+1]
        
        # find expected number of false positives
        fp.result <- binom.false.pos(my.mu, my.alpha, my.beta, Leff, alpha.level=alpha.level, 
                                     by=by, N=N, r=effect.size)
        
        # save binomial data
        tmp.df <- data.frame(sample.size=samp.sizes)
        tmp.df["Power"] <- power.result.bbd
        tmp.df['sample min'] <- bbd.samp.size.min
        tmp.df['sample max'] <- bbd.samp.size.max
        tmp.df['CV'] <- mycv
        tmp.df['alpha.level'] <- alpha.level
        tmp.df['effect.size'] <- effect.size
        tmp.df['mutation.rate'] <- my.mu
        tmp.df["FP"] <- fp.result
        result.df <- rbind(result.df, tmp.df)
      }
      
      # save binomial data
      power.result.binom <- binom.power(my.mu, Leff, alpha.level=alpha.level,
                                        by=by, N=N, r=effect.size)
      binom.samp.size.min <- samp.sizes[min(which(power.result.binom>=desired.power))]
      binom.samp.size.max <- samp.sizes[max(which(power.result.binom<desired.power))+1]
      tmp.df["Power"] <- power.result.binom
      tmp.df['sample min'] <- binom.samp.size.min
      tmp.df['sample max'] <- binom.samp.size.max
      tmp.df['CV'] <- 0
      result.df <- rbind(result.df, tmp.df)
    }
  }
}
# adjust mutation rates back to the average
result.df$mutation.rate <- result.df$mutation.rate / fg
# conver to factor
result.df$mutation.rate <- factor(result.df$mutation.rate, levels=unique(result.df$mutation.rate))
result.df$effect.size <- factor(result.df$effect.size, levels=unique(result.df$effect.size))

###########################
# plot results
###########################
# power if adjust for over-dispersion with beta-binomial model
plot.df <- result.df[result.df['effect.size']==.02,]
plot.df <- result.df[result.df['alpha.level']==1e-4,]
#plot.df <- plot.df[plot.df['mutation.rate']==1.95e-05,]
pwr.plot2 <- ggplot(plot.df, aes(x=sample.size, y=Power, group=CV, color=mutation.rate)) + 
  geom_point() + 
  xlim(0, N) +
  scale_x_continuous(limits=c(0,5000)) +
  facet_wrap(~CV+effect.size) +
  theme(legend.position="top") + 
  geom_hline(aes(yintercept=desired.power))
pwr.plot2 + myTextTheme


# binomial power model
pwr.plot <- ggplot(plot.df, aes(x=sample.size, y=Power, group=effect.size, color=mutation.rate)) + 
  geom_point() + 
  theme(legend.position="top") + 
  xlim(0, N) +
  scale_x_continuous(limits=c(0,5000)) +
  facet_wrap(~effect.size) +
  geom_hline(aes(yintercept=desired.power))
pwr.plot + myTextTheme

# false positives
fp.plot <- ggplot(result.df, aes(x=sample.size, y=FP, group=CV, color=mutation.rate)) + 
  geom_point() + 
  scale_x_continuous(limits=c(0,2000)) +
  scale_y_continuous(limits=c(0,100)) +
  facet_wrap(~CV, ncol=1) + 
  theme(legend.position="top")
fp.plot + myTextTheme

######################
# Save result to text file
######################
file.path <- "/Users/ctokheim/Dropbox (Karchin Lab)/notebooks/prob2020/brainstorming_meeting_6_2015/data/power_false_postives.txt"
write.table(result.df, file.path, sep='\t')

#grid.arrange(pwr.plot, fp.plot, ncol=1)


par(mfrow=c(2,1))
plot(samp.sizes, power.result.bonf, type="l")
lines(samp.sizes, power.result.bh, type="l")
plot(samp.sizes, fp.result.bh, type="l")
lines(samp.sizes, fp.result.bonf, type="l")