# More figures for paper
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

#-------------------------
# (1) Bimodal geyser data
#-------------------------

# The geyser dataset contains the time intervals between eruptions of the Old
# Faithful Geyser

# Figure 1
# Estimate density with the 5 candidate methods
geyser_dat <- multimode::geyser
geyser_res <- estimate_methods(geyser_dat)
plot_methods(geyser_dat, geyser_res, yticks = c(0,0.05))



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

# Figure 2
# Estimate density with the 5 candidate methods
rt_dat <- s01gs$rt
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

# Simulate data from estimated ex-Gaussian
set.seed(79)
exGauss_simdat <- rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1/tau)

# Make GIF to visualize how 5 methods fit the data as more samples become available
#make_gif(exGauss_simdat, "exGauss_simdat", max_y = 0.01, generatingexgauss = c(mu,sigma,tau),
#         xmin = 130, xmax = 530)


# Make plots that show how R density vs Amoroso fit the simulated data at diff n

# Define the sample sizes
n_vec <- c(25, 50, 75, 100)

# Initialize final list
final_list <- list()

# Loop over each sample size
for (n in n_vec) {
  # Extract first n observations
  dat <- exGauss_simdat[1:n]

  # Estimate methods
  res <- estimate_methods(dat)

  # Get maxL amo
  maxL_amo_id <- hellcdf_vs_hellpdf(dat)$maxL_amo

  # Store in a list
  final_list[[paste0("list_", n)]] <- list(
    dat = dat,
    res = res,
    maxL_amo_id = maxL_amo_id
  )
}


##########
# n = 25 #
##########
# Extract object
dat <- final_list$list_25$dat
res <- final_list$list_25$res

# Fit plots
fit_rdens_amo <- plot_rdens_amo(dat, res, xmin = 100, xmax = 600, ymax = 0.012,
                                yticks = c(0,0.012),
                                generatingexgauss = c(mu, sigma, tau))
fit_rdens <- fit_rdens_amo[[1]]
fit_amo <- fit_rdens_amo[[2]]


# Q-Q plots
qq_rdens <- theoretical_qq(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau),
                           rev = F)
qq_amo <- theoretical_qq(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau),
                         rev = F)
# P-P plots
pp_rdens <- theoretical_pp(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau))
pp_amo <- theoretical_pp(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau))
# Arrange them
final_plot_25 <- (fit_rdens + (qq_rdens / pp_rdens) + 
                     fit_amo + (qq_amo / pp_amo)) +
  plot_layout(ncol = 4, widths = c(3, 1, 3, 1)) +
  plot_annotation(
    title = "n = 25",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", size = 20, face = "bold", hjust = 0.5)
    )
  )

# Display final plot
final_plot_25


##########
# n = 50 #
##########
# Extract object
dat <- final_list$list_50$dat
res <- final_list$list_50$res

# Fit plots
fit_rdens_amo <- plot_rdens_amo(dat, res, xmin = 100, xmax = 600, ymax = 0.012,
                                yticks = c(0,0.012),
                                generatingexgauss = c(mu, sigma, tau))
fit_rdens <- fit_rdens_amo[[1]]
fit_amo <- fit_rdens_amo[[2]]

# Q-Q plots
qq_rdens <- theoretical_qq(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau))
qq_amo <- theoretical_qq(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau),
                         rev = TRUE)

# P-P plots
pp_rdens <- theoretical_pp(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau))
pp_amo <- theoretical_pp(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau))

# Arrange them
final_plot_50 <- (fit_rdens + (qq_rdens / pp_rdens) + 
                     fit_amo + (qq_amo / pp_amo)) +
  plot_layout(ncol = 4, widths = c(3, 1, 3, 1)) +
  plot_annotation(
    title = "n = 50",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", size = 20, face = "bold", hjust = 0.5)
    )
  )

# Display final plot
final_plot_50


##########
# n = 75 #
##########
# Extract objects
dat <- final_list$list_75$dat
res <- final_list$list_75$res

# Fit plots
fit_rdens_amo <- plot_rdens_amo(dat, res, xmin = 100, xmax = 600, ymax = 0.012,
                                yticks = c(0,0.012),
                                generatingexgauss = c(mu, sigma, tau))

fit_rdens <- fit_rdens_amo[[1]]
fit_amo <- fit_rdens_amo[[2]]

# Q-Q plots
qq_rdens <- theoretical_qq(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau))
qq_amo <- theoretical_qq(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau),
                         rev = TRUE)

# P-P plots
pp_rdens <- theoretical_pp(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau))
pp_amo <- theoretical_pp(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau))

# Arrange them
final_plot_75 <- (fit_rdens + (qq_rdens / pp_rdens) + 
                     fit_amo + (qq_amo / pp_amo)) +
  plot_layout(ncol = 4, widths = c(3, 1, 3, 1)) +
  plot_annotation(
    title = "n = 75",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", size = 20, face = "bold", hjust = 0.5)
    )
  )

# Display final plot
final_plot_75


###########
# n = 100 #
###########
# Extract object
dat <- final_list$list_100$dat
res <- final_list$list_100$res

# Fit plots
fit_rdens_amo <- plot_rdens_amo(dat, res, xmin = 100, xmax = 600, ymax = 0.012,
                                yticks = c(0,0.012),
                                generatingexgauss = c(mu, sigma, tau))
fit_rdens <- fit_rdens_amo[[1]]
fit_amo <- fit_rdens_amo[[2]]

# Q-Q plots
qq_rdens <- theoretical_qq(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau))
qq_amo <- theoretical_qq(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau),
                         rev = TRUE)

# P-P plots
pp_rdens <- theoretical_pp(dat, res, method_id = "rdens",
                           generatingexgauss = c(mu,sigma,tau))
pp_amo <- theoretical_pp(dat, res, method_id = maxL_amo_id,
                         generatingexgauss = c(mu,sigma,tau))

# Arrange them
final_plot_100 <- (fit_rdens + (qq_rdens / pp_rdens) + 
                     fit_amo + (qq_amo / pp_amo)) +
  plot_layout(ncol = 4, widths = c(3, 1, 3, 1)) +
  plot_annotation(
    title = "n = 100",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", size = 20, face = "bold", hjust = 0.5)
    )
  )

# Display final plot
final_plot_100




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
