################################################################################
### REPRODUCE ERROR IN EX-GAUSSIAN SIMULATIONS
################################################################################

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(R.utils)


# Load functions from Github
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))

# Specify the 3 exGaussian parameter sets
par_exgauss1 <- c(1, 0.1, 2)
par_exgauss2 <- c(1, 0.5, 2)
par_exgauss3 <- c(1, 0.9, 3)

pars_list <- list(par_exgauss1, par_exgauss2, par_exgauss3)


nvec <- c(25,50,100,200)
n <- 50
nrep <- 100
seedvec <- seq(1,nrep)
pars <- pars_list[[1]]
i <- 29
set.seed(seedvec[i])
dat <- rexGAUS(n, mu = pars[1], sigma = pars[2], nu = pars[3])
test <- dat[4]
train <- dat[-4]
dat <- train


safe_execute <- function(expr, object_name, data_vector) {
  env <- new.env()  # Create a new environment
  env$dat <- data_vector  # Assign data_vector to 'dat' in the environment
  
  tryCatch({
    withTimeout({
      eval(bquote(.(expr)), envir = env)  # Use the proper environment
      #eval(bquote(.(quote(scdensity(dat, constraint = "twoInflections+")), "scKDE_2infplus", dat))), envir = env)
    }, timeout = 15)
  }, TimeoutException = function(e) {
    cat(sprintf("Timeout error with fitting %s: Execution took longer than 15 seconds; Other methods were still fit.\n",
                object_name))
    return(NA)  # Assigns NA if execution times out
  }, error = function(e) {
    cat(sprintf("Error with fitting %s: %s; Other methods were still fit.\n",
                object_name, e$message))
    return(NA)  # Assigns NA if another error occurs
  })
}

#-------------------
rdens = safe_execute(quote(density(dat)), "rdens", dat)
scKDE_2infplus = safe_execute(quote(scdensity(dat, constraint = "twoInflections+")), "scKDE_2infplus", dat)
# mnorm = {
#   mnorm_fit <- safe_execute(quote(densityMclust(dat, plot = FALSE)), "mnorm", dat)
#   if (length(mnorm_fit)>1) {
#     mnorm_fit$x <- density(dat)$x  # Use same x as base R density
#     mnorm_fit$y <- predict_mnorm(mnorm_fit$x, mnorm_fit, plot = FALSE)
#   }
#   mnormfit}