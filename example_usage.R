#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#########################################################
### Example usage of the functions in this repository ###
#########################################################

### Source functions
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/dAmoroso.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/plot_amoroso.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso_np.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))



##################
### dAmoroso() ###
##################

# -> get density values of the Amoroso distribution for a set of x values
# and Amoroso parameters

# Define range of x
xvec <- seq(0, 10, length.out = 100)

# Get density values
dAmoroso(xvec, a=4, l=1, c=-7, mu=0)



######################
### plot_amoroso() ###
######################

# -> Plot the Amoroso density for a range of x and a set of 4 parameter values

# Define range of x
xvec <- seq(0, 10, length.out = 1000)

# Plot an Amoroso with default arguments
plot_amoroso(xvec, a=4, l=1, c=-7, mu=0)

# If minimal = FALSE, the plot becomes more informative
plot_amoroso(xvec, a=4, l=1, c=-7, mu=0, minimal=FALSE)

# We can also plot multiple Amorosos in different colours in one plot
plot_amoroso(xvec, a=4, l=1, c=-7, mu=0, minimal=FALSE, col="magenta")
lines(xvec, dgg4(xvec, a=4, l=1, c=-5, mu=0), col="seagreen")
lines(xvec, dgg4(xvec, a=4, l=1, c=-2, mu=0), col="darkorange2")



##########################
### estimate_amoroso() ###
##########################

# -> estimates the Amoroso for a vector of datapoints using the various
# estimation methods descibed in Combes et al.(2022)

# Generate data from normal distribution
set.seed(72)
dat <- rnorm(100, mean = 70, sd = 10)

# Estimate Amorosos on data
res <- estimate_amoroso(dat, plot = TRUE)

# Tibble with all fitted Amoroso models
print(res$all_models)

# Maximum likelihood Amoroso
print(res$max_L_model)

# Minimum BIC Amoroso
print(res$min_BIC_model)

# Generate bimodal data
set.seed(55)
n1 <- rnorm(500, mean = 2, sd = 0.5) # First mode
n2 <- rnorm(500, mean = 7, sd = 1) # Second mode
dat <- c(n1, n2)

# Estimate Amorosos on bimodal data
res <- estimate_amoroso(dat, plot = TRUE)

# -> We see that the Amoroso really struggles with the bimodal dataset
# -> This is because the Amoroso is constrained to be unimodal!



#############################
### estimate_amoroso_np() ###
#############################

# -> Estimate Amoroso and common nonparametric methods for a vector of data
# -> To reduce running time, only 3 Amorosos are estimated:
# MLE, Hellinger PDF and Hellinger CDF
# -> I decided for this compromise (fit only one of the CDF methods and one of
# the PDF methods) because fits from the CDF methods tend to be similar to each
# other and fits from the PDF methods tend to be similar.

# Generate data from normal distribution
set.seed(26)
dat <- rnorm(40, mean = 30, sd = 5)

# Estimate 3 Amorosos and the nonparametric methods
res <- estimate_amoroso_np(dat)

# Returns:
#-> list of all models
glimpse(res$modlist)
#-> list of all valid models (i.e., successfully estimated)
glimpse(res$modlist_valid) 
#-> list of all valid models, interpolated to cover the same x range
glimpse(res$modlist_valid_interp)

# For example, we can extract the x and y values from the Bernstein method
glimpse(res$modlist$bern1$x)
glimpse(res$modlist$bern1$y)


# We can also add the true data-generating normal distribution
res <- estimate_amoroso_np(dat, generatingnormal = c(30,5))

# The methods that fail to fit are skipped
dat <- rexp(100,3)
estimate_amoroso_np(dat)



