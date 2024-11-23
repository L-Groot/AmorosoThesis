# Source basis function for estimating Bernstein from Turnbull & Ghosh (2014)
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/",
              "heads/main/umdz_Rfunction.R"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

estimate_bernstein <- function(dat,
                               breaks = 20,
                               plot = FALSE,
                               n = 512,
                               main = NULL,
                               bound_type = "sd") {
  
  # Extract vector name
  vecname <- paste0("'", deparse(substitute(dat)), "'")
  
  # Remove any missing values
  num_nas_to_remove <- length(dat) - length(na.omit(dat))
  dat <- na.omit(dat)
  
  # Print how many NAs were removed
  if(num_nas_to_remove > 0) {
    message(c("WARNING: ", as.character(num_nas_to_remove), " NAs were removed from the data.\n"))
    cat("--------------------------------------------------------------------\n")
  }
  
  # Estimate Bernstein
  res <- umd(dat, bound.type = bound_type)
  x_min <- res$lower
  x_max <- res$upper
  x_vals <- seq(x_min, x_max, length.out = n)
  y_vals <- res$dumd(x_vals)
  
  # Set xlim for plot a little beyond the range of the x values
  x_span <- range(x_vals)[2]-range(x_vals)[1]
  x_buffer <- 0.05*x_span
  xlim <- c((x_min-x_buffer), (x_max+x_buffer))
  
  # Create main title
  # -> If custom title specified, add that
  if (!is.null(main)) {
    main <- main
  # -> Else add automatic title
  } else {
    main <- paste0("Bernstein Polynomial (bound.type = '", bound_type, "')")
  }
  
  # Optional plot
  if (plot == TRUE) {
    # Plot histogram
    hist(dat, breaks = breaks, xlim = xlim, freq = FALSE, main = main,
         xlab = "Value", ylab = "Density", border = F, col = "#efe5e5")
    # Add Bernstein fit
    lines(x_vals, y_vals, col = "magenta3", lwd = 3)
  }
  
  # Return object
  return(invisible(
      list(res,
           x = x_vals,
           y = y_vals
      )
    )
  )
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Test the function

#set.seed(63)
#dat <- rnorm(1000,2,0.7)

#bern_sd <- estimate_bernstein(dat, bound_type = "sd", plot = TRUE)
#bern_carv <- estimate_bernstein(dat, bound_type = "Carv", plot = TRUE)
