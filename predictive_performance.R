#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Source functions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Test all of the methods for many different kinds of data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

##############################
### Data-generating Normal ###
##############################

### Small n ###
n <- 30
set.seed(29)
dat <- rnorm(n)
hist(dat, breaks = 10)
res <- get_pp(dat)
res$likelihood_tib_avg %>% arrange(desc(logL_avg))



# -> mixed normal very unsmooth at low n
# -> MLE Amoroso regularly fits weird spikes at low n -> CDF and PDF do better




### Predictive performance for normal

n <- 25
mean <- 12
sd <- 3
set.seed(80)
dat <- rnorm(n, mean=100, sd=sd)
hist(dat,breaks=10)
estimate_methods(dat)

n <- 50
sd <- 10
set.seed(80)
dat <- rnorm(n, mean=100, sd=sd)
#hist(dat,breaks=10)
res <- estimate_methods(dat)
#estimate_bernstein(dat, plot = TRUE, bound_type = "Carv")
#estimate_amoroso(dat, plot=T)

n <- 75
sd <- 10
set.seed(80)
dat <- rnorm(n, mean=100, sd=sd)
hist(dat,breaks=10)
estimate_methods(dat)
