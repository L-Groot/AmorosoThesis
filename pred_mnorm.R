# Packages
require(essHist)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Function that predicts density for new data from a mixed normal model
# (coming from densityMclust())

pred_mnorm <- function(x, mod, plot=TRUE, breaks=20) {
  
   # Weight vector (weight of each component)
   propvec <- mod$parameters$pro
   propvec <- propvec/sum(propvec) #normalize to make sure it sums to exactly 1
   # Mean vector
   meanvec <- mod$parameters$mean
   # Variance vector
   varvec <- mod$parameters$variance$sigmasq
   # Model type ("X","E" or "V")
   modeltype <- mod$parameters$variance$modelName
   # If there are n components with equal variance, repeat the variance n times
   if(modeltype == "E") {
     varvec <- rep(varvec, length(meanvec))
   }
   
   # Make predictions
   pred_y <- dmixnorm(x, mean = meanvec, sd = varvec, prob = propvec)
   
   # Optional: Plot
   if (plot == TRUE) {
     # Extract nr of components
     k <- length(meanvec)
     # Extract data
     dat <- mod$data
     # Extract ordered density values
     dens_y <- mod$density[order(dat)]
     # Make plot title with nr of components
     if(k == 1) {
       main <- c("Mixed Normal (k=1)")
     } else {
       main <- paste0("Mixed Normal (k=",k,")")
     }

     # Make histogram
     hist(dat, breaks = breaks, ylim = c(0,max(dens_y)),
          freq = FALSE, col = "lightblue", main = main,
          xlab = "Data", border = "white")
     
     # Add mixed normal fit
     lines(sort(dat), dens_y, col = "red", lwd = 2)
   }
   
   # Return
   return(pred_y)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# TESTS
# set.seed(123)
# 
# # Component 1
# n1 <- 500; mean1 <- 5; sd1 <- 1.2
# 
# # Component 2 (same sd as 1)
# n2 <- 500
# mean2 <- 12
# sd2 <- 1.22
# 
# # Component 3 (diff sd)
# n3 <- 500
# mean3 <- 19
# sd3 <- 4.3
# 
# # Generate the three normal components
# data1 <- rnorm(n1, mean = mean1, sd = sd1)
# data2 <- rnorm(n2, mean = mean2, sd = sd2)
# data3 <- rnorm(n3, mean = mean3, sd = sd3)
# 
# # Create two bimodal datasets
# bimodal_data_equalvar <- c(data1, data2) # -> components have approx. equal var
# bimodal_data_unequalvar <- c(data2, data3) # -> components have unequal var
# 
# # Fit mixed normal
# mod1 <- densityMclust(data1) # -> one component
# mod2 <- densityMclust(bimodal_data_equalvar) # -> two components w. equal var.
# mod3 <- densityMclust(bimodal_data_unequalvar) # -> two components, unequal var.
# 
# # Three new x values
# xvec <- c(3.2,4.5,6.8)
# 
# # Predict densities for three new datapoints
# pred_mnorm(xvec, mod1)
# pred_mnorm(xvec, mod2)
# pred_mnorm(xvec, mod3)

