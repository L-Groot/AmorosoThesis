source("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/heads/main/plot_amoroso.R")
source("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/heads/main/estimate_amoroso.R")

# plot_amoroso()
xx <- seq(-4, 9, by=0.01)
plot_amoroso(xx, a=2, l=0.3, c=5, mu=-3, col="purple")

# estimate_amoroso()
dat <- iris$Sepal.Length
estimate_amoroso(dat, plot = 0)
estimate_amoroso(dat, plot = 1)
estimate_amoroso(dat, plot = 2)
estimate_amoroso(dat, plot = 3)
