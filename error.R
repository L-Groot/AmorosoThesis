################################################################################
### REPRODUCE ERROR IN EX-GAUSSIAN SIMULATIONS
################################################################################

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


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

# Nr of iterations per sample size and distribution
nrep <- 100

# Make vector of seeds
seedvec <- seq(1,nrep)

pars <- pars_list[[1]]

i <- 29

set.seed(seedvec[i])

dat <- rexGAUS(n, mu = pars[1], sigma = pars[2], nu = pars[3])

# Make train and test sets
test <- dat[4]
train <- dat[-4]

estimate_methods(train)

# genpar <- c(pars[1],pars[2],pars[3])
# gendist <- "exgaussian"
# np_methods = c("rdens", "scKDE_2infplus")



# Get predictions for the current test observation
#pp_df <- get_pp_measures(train, test, genpar, gendist, np_methods)
