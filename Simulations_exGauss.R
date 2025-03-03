# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions from Github
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))

# Specify the 3 exGaussian parameter sets
par_exgauss1 <- c(1, 0.1, 2)
par_exgauss2 <- c(1, 0.5, 2)
par_exgauss3 <- c(1, 0.9, 3)

# Put parameter sets in a list
pars_list <- list(par_exgauss1, par_exgauss2, par_exgauss3)

# Vector of sample sizes
nvec <- 25
#nvec <- c(25,50,100)

# Nr of iterations per sample size and distribution
nrep <- 2
#nrep <- 100

# Make vector of seeds
seedvec <- seq(1,nrep)

# Initialize empty list for results
res_exgauss <- list(
  pars1 = list(
    pars = c(),
    n25 = list(),
    n50 = list(),
    n100 = list()
  ),
  pars2 = list(
    pars = c(),
    n25 = list(),
    n50 = list(),
    n100 = list()
  ),
  pars3 = list(
    pars = c(),
    n25 = list(),
    n50 = list(),
    n100 = list()
  )
)

#-------------------------------------------------------------------------------
# Loop through parameter sets
for (parsetnr in 1:length(pars_list)) {
  
  # Get parameters
  pars <- pars_list[[parsetnr]]
  
  # Store parameters in results
  res_exgauss[[parsetnr]]$pars <- pars
  
  # Loop through sample sizes
  for (n_ix in 1:length(nvec)) {
    
    # Get current sample size
    n <- nvec[n_ix]
    
    # Initialize dataframes for results
    # -> Winning models
    win_df <- data.frame(i = 1:nrep,
                         max_logL_method = character(nrep),
                         min_MSE_method = character(nrep))
    # -> NA in medL measure logbook
    na_medL <- data.frame(i = 1:nrep,
                          rdens = numeric(nrep),
                          scKDE_2infplus = numeric(nrep),
                          mnorm = numeric(nrep),
                          amo_hell_cdf = numeric(nrep),
                          amo_hell_pdf = numeric(nrep))
    # -> NA in logL measure logbook
    na_logL <- data.frame(i = 1:nrep,
                          rdens = numeric(nrep),
                          scKDE_2infplus = numeric(nrep),
                          mnorm = numeric(nrep),
                          amo_hell_cdf = numeric(nrep),
                          amo_hell_pdf = numeric(nrep))
    
    
    # Loop through iterations (100 iterations per sample size and parameter set)
    for (i in 1:nrep) {
      
      # Print info of current iteration
      cat("mu = ", pars[1], ", sigma = ", pars[2], ", nu = ", pars[3], " | ",
          "n = ", n, " | ", "i = ", i, "\n")
      
      # Simulate data
      set.seed(seedvec[i])
      dat <- rexGAUS(n, mu = pars[1], sigma = pars[2], nu = pars[3])
      
      # Estimate methods and get PP measures
      res <- get_pp(dat, method = "loocv",
                    generating_exgauss = c(pars[1], pars[2], pars[3]))
      
      # Identify max-L and min-MSE method
      max_logl_meth <- rownames(res)[which.max(ifelse(is.na(res$logL_avg), -Inf, res$logL_avg))]
      min_mse_meth <- rownames(res)[which.min(ifelse(is.na(res$mse_avg), Inf, res$mse_avg))]
      
      # Store them in results dataframe
      win_df$max_logL_method[i] <- max_logl_meth
      win_df$min_MSE_method[i] <- min_mse_meth
      
      # Check if any method has NA for likelihood measures
      methods <- rownames(res)
      na_logL_flags <- as.numeric(is.na(res$logL_avg)) # TRUE if NA, FALSE otherwise
      na_medL_flags <- as.numeric(is.na(res$medL_avg))
      
      # Store the NA flags in results dataframe
      method_match <- intersect(methods, colnames(na_logL))
      na_logL[i, method_match] <- na_logL_flags
      na_medL[i, method_match] <- na_medL_flags
    }
    
    # Add the results dataframes of this n to the results list
    res_exgauss[[parsetnr]][[n_ix+1]] <- list(
      win_df = win_df,
      na_logL = na_logL,
      na_medL = na_medL
    )
    
  }
}

#-------------------------------------------------------------------------------
# Save list!
saveRDS(res_exgauss, file = "res_exgauss.rds")

