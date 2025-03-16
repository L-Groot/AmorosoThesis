par_amo1 <- c(a = 5, l = 0.3, c = 5, mu = 0)
par_exgauss1 <- c(1, 0.1, 2)
par <- par_exgauss1

n <- 2000
set.seed(40)
dat <- rexGAUS(n,mu=par[1],sigma=par[2],nu=par[3])
length(dat)


# Wiggly mixed normal
res_mnorm <- densityMclust(dat, plot = TRUE, G =1999)
df$nr_components[i] <- res_mnorm$G


# Nr of components in mixed normal and wiggly at increasing sample size
increments <- seq(10,2000,by=100)

df <- data.frame(n = increments, nr_components = numeric(20))

for (i in 1:length(increments)) {
  
  x <- increments[i]
  print(x)
  data <- dat[1:x]
  print(length(data))
  res <- estimate_methods(data)
  plot_some_methods(data, res, method_to_plot = "mnorm", lwd = 1.5, legend = F)
  #readline(prompt = "Press [Enter] to continue to the next iteration...")
  res_mnorm <- densityMclust(data, plot = TRUE, G =1000)
  df$nr_components[i] <- res_mnorm$G
  #readline(prompt = "Press [Enter] to continue to the next iteration...")
  
}



res <- densityMclust(dat, G = 2, plot = T)
(nr_comp <- res$G)

print(paste("Iteration:", i))
print(paste("Number of components:", nr_comp))

# Wait for user input before proceeding
#readline(prompt = "Press [Enter] to continue to the next iteration...")


