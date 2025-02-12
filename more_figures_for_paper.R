# More figures for paper
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

#-------------------------
# (1) Bimodal geyser data
#-------------------------

# The geyser dataset contains the time intervals between eruptions of the Old
# Faithful Geyser

# Estimate density with the 5 candidate methods
# geyser_dat <- multimode::geyser
# geyser_res <- estimate_methods(geyser)
# plot_methods(geyser_dat, geyser_res)



#----------------------
# (2) Unimodal RT data
#----------------------

# Get the reaction times from participant 1 from the simple RT task (press key
# as fast as possible upon stimulus presentation)
# Data from "Estimation and Interpretation of 1/f Noise in Human Cognition"
# (Wagenmakers, Farrell, & Ratcliff, 2002).

s01gs <- read.table("noisedat/S01SS.DAT") %>%
  rename(
    key = V1,
    rt = V2,
    ignore = V3,
    odd = V4,
    rsi = V5,
    stim = V6
  ) %>%
  # Compute mean and standard deviation
  mutate(
    rt_mean = mean(rt),
    rt_sd = sd(rt)
  ) %>%
  # Remove RTs outside mean ± 2*SD
  filter(rt >= (rt_mean - 2 * rt_sd) & rt <= (rt_mean + 2 * rt_sd)) %>%
  # Drop temporary columns
  select(-rt_mean, -rt_sd)


# Estimate density with the 5 candidate methods
rt_dat <- s01gs$rt
# hist(rt_dat, breaks = 30)
rt_res <- estimate_methods(rt_dat)
plot_methods(rt_dat, rt_res, yticks = c(0,0.01))

geyser_res <- estimate_methods(geyser)
plot_methods(geyser, geyser_res, yticks = c(0,0.05))



# Estimate ex-gaussian on the data
#install.packages("ExGaussEstim")
library(ExGaussEstim)
exGauss_par <- BayesianExgaussian(length(rt_dat), rt_dat, nSamples = 5000, Ti = 2500)


# Simulate data from ex-gaussian
# Use the estimated parameters
mu <- exGauss_par$mu
sigma <- exGauss_par$sigma
tau <- exGauss_par$tau
# Use same sample size as empirical data
n <- length(rt_dat)

# Simulate data
set.seed(80)
exGauss_simdat <- rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1/tau)
#exGauss_simdat <- exGauss_simdat[1:33]
# Estimate methods on simulated data
#resall <- estimate_methods(exGauss_simdat)
#plot_methods(exGauss_simdat, resall, ymax = 0.012, yticks = c(0,0.012), generatingexgauss = c(mu,sigma,tau))


make_gif(exGauss_simdat, "exGauss_simdat", max_y = 0.01, generatingexgauss = c(mu,sigma,tau),
         xmin = 130, xmax = 530)



#--------------------------------
# (3) Unimodal MCMC samples data
#--------------------------------








#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# bs_and_adjKDE <- function(data, breaks = 20) {
#   
#   par(mfrow=c(1,6))
#   
#   bs_sd <- estimate_bernstein(data, plot=TRUE, breaks = breaks, bound_type = "sd")
#   
#   xlim <- bs_sd$xlim
#   
#   bs_carv <- estimate_bernstein(data, plot=TRUE, breaks = breaks, bound_type = "Carv",
#                                 xlim = xlim)
#   
#   sc <- scdensity(data, constraint = "unimodal")
#   hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "adjKDE-unimodal",
#        xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
#   lines(sc$x, sc$y, col = "magenta3", lwd = 2)
# 
#   sc <- scdensity(data, constraint = "twoInflections")
#   hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "adjKDE-2Inf",
#        xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
#   lines(sc$x, sc$y, col = "magenta3", lwd = 2)
#   
#   sc <- scdensity(data, constraint = "twoInflections+")
#   hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "adjKDE-2Inf+",
#        xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
#   lines(sc$x, sc$y, col = "magenta3", lwd = 2)
#   
#   hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "R density()",
#        xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
#   lines(density(data)$x, density(data)$y, col = "magenta3", lwd = 2)
#   
# }
# 
# # erratic behavior of MLE Amoroso
# # CDF and PDF methods tend to be similar
# set.seed(53)
# data <- rnorm(30)
# estimate_amoroso(data, plot = 2)
# estimate_methods(data, plot = TRUE)
# 
# 
# # undesirable sharp cutoffs from Bernstein sd
# # unreasonably wide fit from Bernstein Carv
# set.seed(28)
# data <- rnorm(50)
# bs_and_adjKDE(data)
