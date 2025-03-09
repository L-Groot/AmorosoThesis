### Figures for presentation

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source functions
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

# Read in empirical RTs from participant 1
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


### All 966 observations
rt_dat <- s01gs$rt

### First 100 observations
rt_dat_100 <- rt_dat[1:100]

# Histogram of first 100
res <- estimate_methods(rt_dat_100)
plot_methods(rt_dat_100, res, method_to_plot = "rdens", alpha = 0, main = "")


estim
estimate_methods(rt_dat_100)



length(rt_dat)
rt_res <- estimate_methods(rt_dat)
plot_methods(rt_dat, rt_res, yticks = c(0,0.01))

# Estimate ex-gaussian on the data
library(ExGaussEstim)
exGauss_par <- BayesianExgaussian(length(rt_dat), rt_dat, nSamples = 5000, Ti = 2500)

# Get the estimated parameters
mu <- exGauss_par$mu
sigma <- exGauss_par$sigma
tau <- exGauss_par$tau

# Use same sample size as original empirical data
n <- length(rt_dat)