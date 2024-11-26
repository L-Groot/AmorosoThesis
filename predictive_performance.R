source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso_np.R"))


### Predictive performance for normal

set.seed(80)
dat <- rnorm(25)
hist(dat,breaks=10)
estimate_amoroso_np(dat)
