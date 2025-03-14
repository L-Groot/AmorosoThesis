par_amo1 <- c(a = 5, l = 0.3, c = 5, mu = 0)
par_exgauss1 <- c(1, 0.1, 2)

n <- 200
dat <- rexGAUS(n,mu=par[1],sigma=par[2],nu=par[3])

par <- par_exgauss1


res <- densityMclust(dat, G = 2, plot = T)
(nr_comp <- res$G)

print(paste("Iteration:", i))
print(paste("Number of components:", nr_comp))

# Wait for user input before proceeding
#readline(prompt = "Press [Enter] to continue to the next iteration...")


