#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Load packages
require(stats)
require(tidyverse)
require(AmoRosoDistrib)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Define function that plots Amoroso for a combination of x range and parameters

plot_amoroso <- function(x_vals, a, l, c, mu,
                         col = "blue", minimal = TRUE,
                         ymax = 1.0,
                         title = TRUE, cex.main = 1.9) {
  
  # Alpha can't be 0
  if (a == 0) {
    stop("alpha cannot be 0")
    
  # If alpha > 0
  } else if (a > 0){
    density_values <- c()
    for (x in x_vals) {
      if(x < mu) { # assign density of 0 to all x smaller than mu
        density_values <- append(density_values,0)
      } else {
        dval <- dgg4(x, a, l, c, mu)
        density_values <- append(density_values, dval)
      }
    }
    
    # If alpha < 0
  } else if (a < 0) {
    density_values <- c()
    for (x in x_vals) {
      if(x > mu) { # assign density of 0 to all x greater than mu
        density_values <- append(density_values,0)
      } else {
        dval <- dgg4(x, a, l, c, mu)
        density_values <- append(density_values, dval)
      }
    }
  }
  
  ####################################
  # Option 1: More informative plots #
  ####################################
  
  if (minimal == FALSE) {
    
    # Plot Amoroso
    plot(x_vals, density_values, type = "l", col = col,
         main = NULL,
         xlab = "", ylab = "",
         lwd = 2,
         ylim = c(0,ymax),
         axes = FALSE)
    
    # Default: Add title with parameters
    if (title == TRUE) {
      title(paste0("Amoroso","(",
                   "α=", a, ", ",
                   "\u2113=", l, ", ",
                   "c=", c, ", ",
                   "μ=", mu,
                   ")",
                  ),
            font = 3, cex.main = cex.main
            )
    }
    
    # Add axes in and labels in more minimal style (for proposal)
    axis(1, cex.axis = 0.8)
    axis(2, cex.axis = 0.8)
    par(las = 0)
    mtext("x", side = 1, line = 2.2, cex = 0.8, font = 1)
    mtext("Density", side = 2, line = 2.6, cex = 0.8, font = 1)
  
    
  ############################################
  # Option 2: Minimalist plot (for proposal) #
  ############################################
    
  } else { # minimal = TRUE
    
    # Plot Amoroso
    plot(x_vals, density_values, type = "l", col = col,
         main = NULL,
         xlab = "", ylab = "",
         lwd = 2,
         ylim = c(0,ymax),
         axes = FALSE)
    
    # Default: Add title with parameters
    if (title == TRUE) {
      title(paste0("α = ", a, ", ",
                   "\u2113 = ", l, ", ",
                   "c = ", c, ", ",
                   "μ = ", mu
                   ),
            font = 3, cex.main = cex.main, line = 0
            )
      axis(1, cex.axis = 1.4)
      par(las = 0)
    }
    
  }
  
}
