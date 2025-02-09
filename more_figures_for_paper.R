# More figures for paper
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

#-------------------------
# (1) Bimodal geyser data
#-------------------------

geyser_dat <- multimode::geyser
geyser_res <- estimate_methods(geyser)
plot_methods(geyser_dat, geyser_res)



#----------------------
# (2) Unimodal RT data
#----------------------

# Data from subject number 1 on the 
s01gs <- read.table("noisedat/S01GS.DAT") %>%
  # Rename columns
  rename(
    key = V1,
    rt = V2,
    ignore = V3,
    odd = V4,
    rsi = V5,
    stim = V6
  ) %>%
  # Remove rows with rt < 500ms or > 2500ms
  filter(rt >= 500 & rt <= 2500)


rt_dat <- s01gs
rt_res <- estimate_methods(rt_dat)


#--------------------------------
# (3) Unimodal MCMC samples data
#--------------------------------







bs_and_adjKDE <- function(data, breaks = 20) {
  
  par(mfrow=c(1,6))
  
  bs_sd <- estimate_bernstein(data, plot=TRUE, breaks = breaks, bound_type = "sd")
  
  xlim <- bs_sd$xlim
  
  bs_carv <- estimate_bernstein(data, plot=TRUE, breaks = breaks, bound_type = "Carv",
                                xlim = xlim)
  
  sc <- scdensity(data, constraint = "unimodal")
  hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "adjKDE-unimodal",
       xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
  lines(sc$x, sc$y, col = "magenta3", lwd = 2)

  sc <- scdensity(data, constraint = "twoInflections")
  hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "adjKDE-2Inf",
       xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
  lines(sc$x, sc$y, col = "magenta3", lwd = 2)
  
  sc <- scdensity(data, constraint = "twoInflections+")
  hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "adjKDE-2Inf+",
       xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
  lines(sc$x, sc$y, col = "magenta3", lwd = 2)
  
  hist(data, breaks = breaks, xlim = xlim, freq = FALSE, main = "R density()",
       xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
  lines(density(data)$x, density(data)$y, col = "magenta3", lwd = 2)
  
}

# erratic behavior of MLE Amoroso
# CDF and PDF methods tend to be similar
set.seed(53)
data <- rnorm(30)
estimate_amoroso(data, plot = 2)
estimate_methods(data, plot = TRUE)


# undesirable sharp cutoffs from Bernstein sd
# unreasonably wide fit from Bernstein Carv
set.seed(28)
data <- rnorm(50)
bs_and_adjKDE(data)
