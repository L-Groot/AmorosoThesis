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


### All 972 observations
rt_dat <- s01gs$rt

### First 100 observations
rt_dat_100 <- rt_dat[1:100]

# Histogram of first 100
res <- estimate_methods(rt_dat_100)
plot_methods(rt_dat_100, res, method_to_plot = "rdens", alpha = 0, main = "",
             xmin = 150, xmax = 550, ymax = 0.017, yticks = c(0,0.017))

# Density fit for first 100
plot_methods(rt_dat_100, res, method_to_plot = "rdens", alpha = 1, main = "",
             xmin = 150, xmax = 550, ymax = 0.017, yticks = c(0,0.017))

# First 250, 500 and all
rt_dat_250 <- rt_dat[1:250]
res <- estimate_methods(rt_dat_250)
plot_methods(rt_dat_250, res, method_to_plot = "rdens", alpha = 1, main = "",
             xmin = 150, xmax = 550, ymax = 0.017, yticks = c(0,0.017))


rt_dat_500 <- rt_dat[1:500]
res <- estimate_methods(rt_dat_500)
plot_methods(rt_dat_500, res, method_to_plot = "rdens", alpha = 1, main = "",
             xmin = 150, xmax = 550, ymax = 0.017, yticks = c(0,0.017))


res <- estimate_methods(rt_dat)
plot_methods(rt_dat, res, method_to_plot = "rdens", alpha = 1, main = "",
             xmin = 150, xmax = 550, ymax = 0.017, yticks = c(0,0.017))

# Estimate ex-gaussian on the data
library(ExGaussEstim)
exGauss_par <- BayesianExgaussian(length(rt_dat), rt_dat, nSamples = 5000, Ti = 2500)

# Get the estimated parameters
mu <- exGauss_par$mu
sigma <- exGauss_par$sigma
tau <- exGauss_par$tau

# Use same sample size as original empirical data
n <- length(rt_dat)

# Simulate new dataset of same size
set.seed(17)
exGauss_simdat <- rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1/tau)

# First 100
simdat_100 <- exGauss_simdat[1:100]
res <- estimate_methods(simdat_100)

# Data-generating distribution
plot_methods(simdat_100, res, method_to_plot = "rdens", generatingexgauss = c(mu,sigma,tau),
             alpha = 0, histfill = "grey100", histoutline = "grey100", main = "",
             xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))

# Histogram of data
plot_methods(simdat_100, res, method_to_plot = "rdens", generatingexgauss = c(mu,sigma,tau),
             alpha = 0, histfill = "grey90", histoutline = "grey80", main = "",
             xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))

# Method fits to first 100
plot_methods(simdat_100, res, method_to_plot = "rdens", generatingexgauss = c(mu,sigma,tau),
             alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
             xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))
plot_methods(simdat_100, res, method_to_plot = "scKDE_2infplus", generatingexgauss = c(mu,sigma,tau),
             alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
             xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))
plot_methods(simdat_100, res, method_to_plot = "amo_hell_cdf", generatingexgauss = c(mu,sigma,tau),
             alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
             xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))
plot_methods(simdat_100, res, method_to_plot = "amo_hell_pdf", generatingexgauss = c(mu,sigma,tau),
             alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
             xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))
plot_methods(simdat_100, res, method_to_plot = "mnorm", generatingexgauss = c(mu,sigma,tau),
             alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
             xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))

# Repeat for 10 other randomly generated 100 observations 
# Create the directory if it doesn't exist
if (!dir.exists("figures1/first_10_sets")) {
  dir.create("figures1/first_10_sets")
}


# Not cherrypicking! Show density() fits for first 10 simulated datasets
seeds <- seq(1, 10)

for (seed in seeds) {
  n <- 100
  set.seed(seed)  # Ensure the random seed is set
  exGauss_simdat <- rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1/tau)
  res <- estimate_methods(exGauss_simdat)
  
  # Create the plot
  p <- plot_methods(exGauss_simdat, res, method_to_plot = "rdens", generatingexgauss = c(mu, sigma, tau),
                    alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
                    xmin = 100, xmax = 600, ymax = 0.011, yticks = c(0,0.011))
  
  # Save the plot to a file
  plot_filename <- paste0("figures1/first_10_sets/set", seed, ".png")  # Save as .png (you can change the format)
  ggsave(plot_filename, plot = p, width = 6, height = 6)  # Adjust the width and height as needed
}


# Explain likelihood by comparing Amoroso vs R density()
plot_some_methods(simdat_100, res, method_to_plot = c("rdens","amo_hell_cdf"), generatingexgauss = c(mu,sigma,tau),
                  alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
                  xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))


plot_some_methods(simdat_100, res, method_to_plot = c("rdens","amo_hell_cdf"),
                  alpha = 1, histfill = "grey90", histoutline = "grey80", main = "",
                  xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))

plot_some_methods(simdat_100, res, method_to_plot = c("amo_hell_cdf"),
                  alpha = 1, histfill = "grey100", histoutline = "grey100", main = "",
                  xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))


plot_some_methods(simdat_100, res, method_to_plot = c("rdens","amo_hell_cdf"),
                  alpha = 0, histfill = "grey90", histoutline = "grey80", main = "",
                  xmin = 100, xmax = 600, ymax = 0.01, yticks = c(0,0.01))






  




