par_amo1 <- c(a = 5, l = 0.3, c = 5, mu = 0)
par_exgauss1 <- c(1, 0.1, 2)

n <- 200
seeds <- 1:20

par <- par_exgauss1

for (i in 1:20) {
  set.seed(seeds[i])
  #dat <- rgg4(n, a = par[1], l = par[2], c = par[3], mu = par[4], 
  #            sequence = FALSE)
  dat <- rexGAUS(n,mu=par[1],sigma=par[2],nu=par[3])
  
  res <- densityMclust(dat, G = 90)
  nr_comp <- res$G
  
  print(paste("Iteration:", i))
  print(paste("Number of components:", nr_comp))
  
  # Wait for user input before proceeding
  readline(prompt = "Press [Enter] to continue to the next iteration...")
}

