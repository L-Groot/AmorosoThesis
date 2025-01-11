### Trying out the BayesMultiMode package

#install.packages("BayesMultiMode")
library(BayesMultiMode)

# Load data
#y <- galaxy
#y <- BayesMultiMode::cyclone$max_wind

## Unimodal (data-gen normal)
# -> density() fit is far from that (wiggly and spurious mode(s))
#set.seed(25)
y <- rnorm(100)

set.seed(25)
y <- rnorm(100)

## Bimodal
#set.seed(123)
# data1 <- rnorm(100, mean = 5, sd = 1)
# data2 <- rnorm(100, mean = 10, sd = 1)
# y <- c(data1, data2)

## Exponential
# -> Amoroso also looks promising here
# lambda <- 0.5 #rate
# n <- 100
# y <- rexp(n, rate = lambda)

## Normal
y <- rnorm(100)


# Estimate univariate mixture with unknown nr of components
set.seed(123) #seed for reproducibility

mix_mcmc <- bayes_fit(data = y,
                             K = 10,
                             dist = "normal",
                             nb_iter = 2000, #default
                             burnin = 1000) #default

# Get mcmc draws (excl burnin)
mix_mcmc$mcmc

summary(mix_mcmc)
# Plot estimated density for 100s draws
plot(mix_mcmc)

# Detect modes with fixed-point algorithm (since we estimated Gaussian mixture)
mode_mcmc <- bayes_mode(mix_mcmc)
summary(mode_mcmc)

# Get probability that density is unimodal
(p_unimodal <- mode_mcmc$p1) #probability of the density being unimodal

# Plot various NP and P methods
# res <- estimate_methods(y)

#broom.mixed
#tidy()
