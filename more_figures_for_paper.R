# More figures for paper

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
