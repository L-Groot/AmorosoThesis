# More figures for paper # Yay, Visualisations!
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/make_gif.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/mnorm_functions.R"))

library(patchwork)

#-------------------------
# (1) Bimodal geyser data
#-------------------------

# The geyser dataset contains the time intervals between eruptions of the Old
# Faithful Geyser

# Figure 1
# Estimate density with the 5 candidate methods
geyser_dat <- multimode::geyser
geyser_res <- estimate_methods(geyser_dat)
plot_methods(geyser_dat, geyser_res,
             xticks = c(30,70,110), yticks = c(0,0.05),
             xmin = 30, xmax = 110,
             cex_main = 18, cex_axislab = 15, cex_axistick = 15)

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
plot_methods(rt_dat, rt_res,
             xticks = c(150,350,550), yticks = c(0,0.01),
             xmin = 150, xmax = 550,
             cex_main = 18, cex_axislab = 15, cex_axistick = 15)

# Estimate ex-gaussian on the data
library(ExGaussEstim)
exGauss_par <- BayesianExgaussian(length(rt_dat), rt_dat, nSamples = 5000, Ti = 2500)

# Get the estimated parameters
mu <- exGauss_par$mu
sigma <- exGauss_par$sigma
tau <- exGauss_par$tau

# Use same sample size as original empirical data
n <- 2000

# Simulate data from estimated ex-Gaussian
set.seed(80)
exGauss_simdat <- rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1/tau)

# Make GIF to visualize how 5 methods fit the data as more samples become available
make_gif(exGauss_simdat, "gif_exgauss_1", max_y = 0.01,
         generatingexgauss = c(mu,sigma,tau),
         xmin = 130, xmax = 530,
         batch_size = 25,
         main_datagen = "Data-generating ex-Gaussian(292,9, 44.3, 37.0), n = 2000")

# Make GIF for a more skewed ex-Gaussian
mu <- 1
sigma <- 0.1
tau <- 2

n <- 2000

set.seed(80)
exGauss_simdat <- rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1/tau)
range(exGauss_simdat)

make_gif(exGauss_simdat, "gif_exgauss_2", max_y = 0.5,
         generatingexgauss = c(mu,sigma,tau),
         xmin = 0, xmax = 15,
         batch_size = 25,
         main_datagen = "Data-generating ex-Gaussian(1, 0.1, 2), n = 2000")

# Make GIF to visualize CDF and PDF Amoroso can approximate the normal
n <- 200
par_norm <- c(100,15)
par <- par_norm

set.seed(40)
dat <- rnorm(n,par[1],par[2])

make_gif(dat, "gif_norm_1", max_y = 0.03,
         generatingnormal = c(100,15),
         xmin = 30, xmax = 180,
         main_datagen = "Data-generating Normal(100, 15), n = 200"
)

# Make GIF for mixed normal
n <- 2000
par_mnorm <- c(0.5,3.7,0.7,6,0.8)
par <- par_mnorm

set.seed(22)
dat <- rmixnorm(n, par[1], par[2], par[3], par[4], par[5])
hist(dat)

#res <- estimate_methods(dat)
#plot_methods(dat,res)

make_gif(dat, "gif_mnorm_1", max_y = 0.3,
         generatingmnorm = c(par[1],par[2],par[3],par[4],par[5]),
         xmin = 0, xmax = 10,
         batch_size = 25,
         main_datagen = "Data-generating Mixed Normal(p1=0.5, mu1=3.7, sd1=0.7, mu2=6, sd2=0.8), n = 2000"
)


# Make GIF 2 for mixed normal
n <- 2000
par_mnorm1 <- c(0.6, 3, 1.2, 7, 1)
par <- par_mnorm

set.seed(22)
dat <- rmixnorm(n, par[1], par[2], par[3], par[4], par[5])
hist(dat)

make_gif(dat, "gif_mnorm_2", max_y = 0.3,
         generatingmnorm = c(par[1],par[2],par[3],par[4],par[5]),
         xmin = 0, xmax = 10,
         batch_size = 25,
         main_datagen = "Data-generating Mixed Normal(p1=0.6, mu1=3, sd1=1.2, mu2=7, sd2=1), n = 2000"
)


# Make GIF for Amoroso
n <- 2000
par_amo1 <- c(a = 2, l = 0.3, c = 5, mu = 0)
par <- par_amo1
set.seed(230)
dat <- rgg4(n, par[1], par[2], par[3], par[4])
hist(dat, freq = F)

res <- estimate_methods(dat)
plot_methods(dat,res)

make_gif(dat, "gif_amo_1", max_y = 0.6,
         generatingamoroso = c(par[1],par[2],par[3],par[4]),
         xmin = 0, xmax = 4,
         batch_size = 25,
         main_datagen = "Data-generating Amoroso(a=2, l=0.3, c=5, mu=0), n = 2000"
)



# Make plots that show how R density vs Amoroso fit the simulated data at diff n

# Define the sample sizes
n_vec <- c(25, 50, 75, 100)

# Initialize final list
final_list <- list()

# Loop over each sample size
# for (n in n_vec) {
#   # Extract first n observations
#   dat <- exGauss_simdat[1:n]
# 
#   # Estimate methods
#   res <- estimate_methods(dat)
# 
#   # Get maxL amo
#   maxL_amo_id <- hellcdf_vs_hellpdf(dat)$maxL_amo
# 
#   # Store in a list
#   final_list[[paste0("list_", n)]] <- list(
#     dat = dat,
#     res = res,
#     maxL_amo_id = maxL_amo_id
#   )
# }


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


#-------------------------------------------------
# Appendix D: CDF and PDF Amorosos look alike
#-------------------------------------------------
set.seed(77)
dat <- rnorm(100, plot = 1)
estimate_amoroso(dat)

#-------------------------------------------------
# Appendix F: Number of Components in Mixed Normal
#-------------------------------------------------
par_exgauss1 <- c(1, 0.1, 2)
par <- par_exgauss1

n <- 2000
set.seed(40)
dat <- rexGAUS(n,mu=par[1],sigma=par[2],nu=par[3])
length(dat)

res <- estimate_methods(dat, G=1910)

plot_some_methods(dat, res,
                  generatingexgauss = c(1,0.1,2),
                  method_to_plot = "mnorm",
                  lwd=0.5,
                  legend = FALSE,
                  bins = 60,
                  cex_main = 17,
                  cex_axislab = 14,
                  cex_axistick = 14,
                  ymax = 0.9,
                  xmin = 0,
                  xmax=20,
                  xticks = c(0,5,10,15,20),
                  yticks = c(0,0.8),
                  main = expression(bold("Mixed Normal with "*bolditalic(G)*"=1910 Components")))


# Nr of components in mixed normal and wiggly at increasing sample size
increments <- seq(10,2000,by=100)

df <- data.frame(n = increments, nr_components = numeric(20))

for (i in 1:length(increments)) {
  x <- increments[i]
  print(x)
  data <- dat[1:x]
  print(length(data))
  #res <- estimate_methods(data)
  #plot_some_methods(data, res, method_to_plot = "mnorm", lwd = 1.5, legend = F)
  #readline(prompt = "Press [Enter] to continue to the next iteration...")
  res_mnorm <- densityMclust(data, plot = TRUE)
  df$nr_components[i] <- res_mnorm$G
  #readline(prompt = "Press [Enter] to continue to the next iteration...")
}

print(df)






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


list.files()
