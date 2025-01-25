source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))

# Convenience function to make histogram
make_hist <- function() {
  
  hist(x, xlim = c(xmin, xmax), ylim = c(0,ymax),
       probability = TRUE, main = "", breaks = 20, axes = FALSE,
       col = "grey95", border = "grey85")
  axis(1)
  axis(2, las = 1)
  
}

# Simulate data from standard normal
n <- 100
mean = 175
sd <- 7
set.seed(12)
x <- rnorm(n, mean, sd)
#x <- x[1:80]

# Fit methods
rdens <- density(x)
amo <- estimate_amoroso(x, plot = F)
mnorm <- densityMclust(x, plot=F)
adjKDE <- scdensity(x, constraint = "twoInflections+")

#################
### Plot fits ###
#################

# Set x range for plots
xmin <- 150
xmax <- 200
#xx <- seq(min(density(x)$x), max(density(x)$x), length = n)
xx <- seq(xmin, xmax, length = 512)

# Set y range for plots based on R density() fit
ybuffer <- 0.6*(range(rdens$y)[2]-range(rdens$y)[1])
ymax <- max(rdens$y)+ybuffer

# Settings
# -> single plot
par(mfrow = c(1,1), cex.main = 1.4, mar = c(4,4,0,0), cex.axis = 1,
    cex.lab = 1.2, bty = "n", font.lab = 2,
    family = "Times New Roman")
# -> 2x3 grid
# par(mfrow = c(2,3), oma = c(2, 2, 2, 2), mar = c(4,4,1,1), cex.axis = 1,
#     cex.main = 1.4, cex.lab = 1.2, bty = "n", font.lab = 2,
#     family = "Times New Roman")

# Plot R density fit
make_hist()
lines(xx, dnorm(xx, mean, sd), col = "grey", lty = 2, lwd = 2)
lines(rdens$x, rdens$y, col = "black", lty = 1, lwd = 2)
#lines(xx, dnorm(xx, mean, sd), col = "grey", lty = 2, lwd = 2)
#points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]),
#       type = "l", lwd = 2, col = "magenta2", lty = 1)

# Plot mixed normal
xy_ordered <- data.frame(x=mnorm$data,y=mnorm$density) %>% arrange(x)
mnorm$x <- xy_ordered$x
mnorm$y <- xy_ordered$y

make_hist()
lines(xx, dnorm(xx, mean, sd), col = "grey", lty = 2, lwd = 2)
points(mnorm$x, mnorm$y,
       type = "l", lwd = 2, col = "black", lty = 1)

# Plot adjusted KDE
make_hist()
lines(xx, dnorm(xx, mean, sd), col = "grey", lty = 2, lwd = 2)
points(adjKDE$x, adjKDE$y,
       type = "l", lwd = 2, col = "black", lty = 1)

# Plot Amoroso MLE fit
pars <- amo$min_BIC_models %>%
  filter(method_ID == "MLE") %>%
  select(a,l,c,mu) %>%
  as.numeric()

make_hist()
lines(xx, dnorm(xx, mean, sd), col = "grey", lty = 2, lwd = 2)
points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]),
       type = "l", lwd = 2, col = "black", lty = 1)

# Plot Amoroso Hell-CDF fit
pars <- amo$min_BIC_models %>%
  filter(method_ID == "HELL-CDF") %>%
  select(a,l,c,mu) %>%
  as.numeric()

make_hist()
lines(xx, dnorm(xx, mean, sd), col = "grey", lty = 2, lwd = 2)
points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]),
       type = "l", lwd = 2, col = "black", lty = 1)

# Plot Amoroso Hell-PDF fit
pars <- amo$min_BIC_models %>%
  filter(method_ID == "HELL-PDF") %>%
  select(a,l,c,mu) %>%
  as.numeric()

make_hist()
lines(xx, dnorm(xx, mean, sd), col = "grey", lty = 2, lwd = 2)
points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]),
       type = "l", lwd = 2, col = "black", lty = 1)
