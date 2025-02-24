#-------------------------------------------------------------------------------
# Mixed normal PDF
mixed_normal_pdf <- function(x, propvec, meanvec, varvec) {
  pdf_vals <- numeric(length(x))
  for (i in seq_along(propvec)) {
    pdf_vals <- pdf_vals + propvec[i] * dnorm(x, mean = meanvec[i], sd = sqrt(varvec[i]))
  }
  return(pdf_vals)
}

#-------------------------------------------------------------------------------
# Mixed normal CDF (for two components)
dmixnorm <- function(x, p, mu1, sd1, mu2, sd2) {
  p * dnorm(x, mean = mu1, sd = sd1) + (1 - p) * dnorm(x, mean = mu2, sd = sd2)
}
# -> p is the proportion (weight) of first component 

#-------------------------------------------------------------------------------
# Function that predicts new data from a densityMclust() mixed normal model

predict_mnorm <- function(x, mod, plot=TRUE, breaks=20) {

  # Weight vector (weight of each component)
  propvec <- mod$parameters$pro
  propvec <- propvec / sum(propvec) # Normalize to make sure it sums to exactly 1
  # Mean vector
  meanvec <- mod$parameters$mean
  # Variance vector
  varvec <- mod$parameters$variance$sigmasq
  # Model type ("X","E" or "V")
  modeltype <- mod$parameters$variance$modelName
  # If there are n components with equal variance, repeat the variance n times
  if (modeltype == "E") {
    varvec <- rep(varvec, length(meanvec))
  }
  
  # Make predictions for x values
  pred_y <- mixed_normal_pdf(x, propvec, meanvec, varvec)
  
  # Optional: Plot
  if (plot == TRUE) {
    k <- length(meanvec) # nr components
    dat <- mod$data
    
    # Calculate x values for plotting the PDF
    x_vals <- seq(min(dat), max(dat), length.out = 1000)
    pdf_vals <- mixed_normal_pdf(x_vals, propvec, meanvec, varvec)
    
    # Make plot title with number of components
    if (k == 1) {
      main <- "Mixed Normal (k=1)"
    } else {
      main <- paste0("Mixed Normal (k=", k, ")")
    }
    
    # Make histogram
    hist(dat, breaks = 20, ylim = c(0, max(pdf_vals, hist(dat, breaks = 20, plot = FALSE)$density)),
         freq = FALSE, col = "lightblue", main = main,
         xlab = "Data", border = "white")
    
    # Add mixed normal PDF
    lines(x_vals, pdf_vals, col = "red", lwd = 2)
  }
  
  # Return
  return(pred_y)
}

