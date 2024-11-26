
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso.R"))

source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_bernstein.R"))

source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso_np.R"))


### Predictive performance for normal

n <- 25
sd <- 10
set.seed(80)
dat <- rnorm(n, mean=100, sd=sd)
hist(dat,breaks=10)
estimate_amoroso_np(dat)

n <- 50
sd <- 10
set.seed(80)
dat <- rnorm(n, mean=100, sd=sd)
#hist(dat,breaks=10)
estimate_amoroso_np(dat)
#estimate_bernstein(dat, plot = TRUE)
#estimate_amoroso(dat, plot=T)

n <- 75
sd <- 10
set.seed(80)
dat <- rnorm(n, mean=100, sd=sd)
hist(dat,breaks=10)
estimate_amoroso_np(dat)
