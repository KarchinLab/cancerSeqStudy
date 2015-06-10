# just check out curves
library(VGAM)
library(reshape2)
library(ggplot2)
library(gridExtra)
pi <- c(.5e-6, 1e-6, 3e-6, 5e-6, 10e-6)   # 2.85e-6  # used genome wide val from a run of mutsig
nonsilentFactor <- 3/4  # assume only 3/4 of mutations are non-silent like in paper
#mu = 3e-6
cv <- .05  # coefficient of variation
ab <- pi * (1-pi) / (cv*pi)^2 - 1
alpha <- pi * ab
beta <- (1-pi)*ab
L <- 1500  # same length as used in paper
Leff <- L * nonsilentFactor
numGenes <- 18500  # reasonable number of genes to assume are being tested
numSamples <- 8251  # number of samples in my data
basesAtRisk <- Leff*numSamples
z <- 1:300


myTextTheme <- theme(axis.text.x = element_text(size=16, angle=-90),
                     axis.text.y = element_text(size=16),  
                     title = element_text(size=24),
                     axis.title.x = element_text(size=20),
                     axis.title.y = element_text(size=20))  # change text size

# plot pmf of beta-binomial
mb <- 1e6
df <- data.frame(mutation.count=z)
col.names <- paste(pi*mb, "/ Mb", sep=" ")
for (i in 1:length(pi)){ 
  pmf <- dbetabinom.ab(z, basesAtRisk, alpha[i], beta[i])
  df[col.names[i]] <- pmf
}
mybbd.dist <- melt(df, 
                   id.vars=c("mutation.count"), 
                   value.name="PMF", 
                   variable.name="mutation.rate")
bbd.plot <- ggplot(data=mybbd.dist, 
                   aes(x=mutation.count, y=PMF, 
                       group=mutation.rate, color=mutation.rate)) + 
            geom_line(size=1.5) + 
            ggtitle("Beta-Binomial (over-dispersion)") + 
            myTextTheme

# binomial distribution distribution
mb <- 1e6
tmp <- data.frame(mutation.count=z)
col.names <- paste(pi*mb, "/ Mb", sep=" ")
for (i in 1:length(pi)){ 
  pmf <- dbinom(z, basesAtRisk, pi[i])
  tmp[col.names[i]] <- pmf
}
mybinom.dist <- melt(tmp, 
                     id.vars=c("mutation.count"), 
                     value.name="PMF", 
                     variable.name="mutation.rate")
binom.plot <- ggplot(data=mybinom.dist, 
                     aes(x=mutation.count, y=PMF, 
                         group=mutation.rate, color=mutation.rate)) + 
              geom_line(size=1.5) + 
              ggtitle("Binomial") + 
              myTextTheme

# plot beta-binomial and binomial PMF's
grid.arrange(bbd.plot, binom.plot, nrow=2)

# plot beta prior for each of the mutation rates
r <- seq(1e-7, 5e-5, by=1e-8)
tmp <- data.frame(mutation.rate=r)
col.names <- paste(pi*mb, "/ Mb", sep=" ")
for (i in 1:length(pi)){ 
  pmf <- dbeta(r, alpha[i], beta[i])
  tmp[col.names[i]] <- pmf
}
beta.dist <- melt(tmp, 
                  id.vars=c("mutation.rate"), 
                  value.name="Density", 
                  variable.name="mean.mutation.rate")
beta.plot <- ggplot(data=beta.dist, 
                    aes(x=log10(mutation.rate), y=Density, 
                        group=mean.mutation.rate, color=mean.mutation.rate)) + 
             geom_line(size=1.5) + 
             ggtitle(paste("Variability in Mutation Rate if CV=", cv)) + 
             myTextTheme
beta.plot
