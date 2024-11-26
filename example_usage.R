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
plot_amoroso(xvec, a=4, l=1, c=-7, mu=0, minimal=F, col="magenta", title = "")
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

# If one method has a spike the y axis is adjusted to cut that off in the plot
dat <- rexp(100,3)
estimate_amoroso_np(dat)

# Any estimation methods that fail to fit are skipped
set.seed(70)
dat <- rgg4(50,4,1,-0.1,0) # produce weird data
hist(dat) # Inspect it
estimate_amoroso_np(dat) # Try to fit Amorosos and NP methods



################
### get_pp() ###
################

# -> Assess predictive performance of the 3 Amoroso and various nonparametric fits
# to a certain dataset

# Generate data from normal distribution
set.seed(93)
dat <- rnorm(50, mean = 30, sd = 7)

# Get predictive performance of Amoroso and NP Methods
res <- get_pp(dat)
# -> default uses k-fold CV with k=5

# The results contain two likelihood measures:
res$likelihood_tib_avg
# -> each model's average log-likelihood of the test data across the folds
# -> each model's average median likelihood of the test data across the folds

# We can also vary k
res <- get_pp(dat, k = 4)

# Instead of k-fold, we can also just do a 'one-shot' train-test split
res <- get_pp(dat, method = "split-half", prop_train = 0.7)
# -> 70% of data is used for training, 30% for testing

# The results are now simply the two likelihood measures on the test set
res$likelihood_tib

# If the data-generating distribution is an Amoroso with known parameters,
# we can also express the likelihood measures as a proportion of the 'true'
# likelihood of the test set.
# This makes it easier to gauge the size of the differences in predictive
# performance between the methods
set.seed(65)
dat <- rgg4(50, a=4, l=1, c=-5, mu=0) # Generate data from Amoroso
res <- get_pp(dat, generating_amoroso = c(4,1,-5,0))
res$likelihood_tib_avg




 

