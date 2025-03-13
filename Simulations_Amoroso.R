################################################################################
### Simulations: Data-generating Amoroso with Restart Capability
################################################################################

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions from Github
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))

# Define parameters for three Amoroso distributions
par_amo1 <- c(a = 2, l = 0.3, c = 5, mu = 0)
par_amo2 <- c(a = 5, l = 0.3, c = 5, mu = 0)
par_amo3 <- c(a = 7, l = 1.3, c = 8, mu = 0)

# Put parameter sets in a list
pars_list <- list(par_amo1, par_amo2, par_amo3)

# Vector of sample sizes
nvec <- c(25, 50, 100, 200)

# Number of iterations per sample size and distribution
nrep <- 100

# Make vector of seeds
seedvec <- seq(1, nrep)

# Specify where to restart
restart_parsetnr <- 1  # The parameter set number to resume from
restart_n <- 25       # The sample size to resume from
restart_i <- 1       # The iteration number to resume from

# Load previous results if available
if (file.exists("res_amo.rds")) {
  res_amo <- readRDS("res_amo.rds")
} else {
  # Initialize empty list for results
  res_amo <- list(
    pars1 = list(pars = c(), n25 = list(), n50 = list(), n100 = list(), n200 = list()),
    pars2 = list(pars = c(), n25 = list(), n50 = list(), n100 = list(), n200 = list()),
    pars3 = list(pars = c(), n25 = list(), n50 = list(), n100 = list(), n200 = list())
  )
}



#-------------------------------------------------------------------------------
# Loop through parameter sets, skipping completed ones
for (parsetnr in restart_parsetnr:length(pars_list)) {
  
  # Get parameters
  pars <- pars_list[[parsetnr]]
  
  # Store parameters in results
  res_amo[[parsetnr]]$pars <- pars
  
  # Loop through sample sizes, skipping completed ones
  for (n_ix in which(nvec == restart_n):length(nvec)) {
    
    # Get current sample size
    n <- nvec[n_ix]
    
    # Initialize dataframes for results if not already present
    if (!("win_df" %in% names(res_amo[[parsetnr]][[n_ix+1]]))) {
      win_df <- data.frame(i = 1:nrep,
                           max_logL_method = character(nrep),
                           min_MSE_method = character(nrep))
      na_medL <- data.frame(i = 1:nrep,
                            rdens = numeric(nrep),
                            scKDE_2infplus = numeric(nrep),
                            mnorm = numeric(nrep),
                            amo_hell_cdf = numeric(nrep),
                            amo_hell_pdf = numeric(nrep))
      na_logL <- data.frame(i = 1:nrep,
                            scKDE_2infplus = numeric(nrep),
                            mnorm = numeric(nrep),
                            amo_hell_cdf = numeric(nrep),
                            amo_hell_pdf = numeric(nrep))
    } else {
      # Load existing results if they exist
      win_df <- res_amo[[parsetnr]][[n_ix+1]]$win_df
      na_medL <- res_amo[[parsetnr]][[n_ix+1]]$na_medL
      na_logL <- res_amo[[parsetnr]][[n_ix+1]]$na_logL
    }
    
    # Loop through iterations, skipping completed ones
    for (i in restart_i:nrep) {
      
      # Print info of current iteration
      cat("Restarting at: a = ", pars[1], ", l = ", pars[2], ", c = ", pars[3],
          ", mu = ", pars[4], " | ", "n = ", n, " | ", "i = ", i, "\n")
      
      # Simulate data
      set.seed(seedvec[i])
      dat <- rgg4(n, a = pars[1], l = pars[2], c = pars[3], mu = pars[4], sequence = FALSE)
      
      # Estimate methods and get PP measures
      res <- get_pp(dat, method = "loocv",
                    generating_amoroso = c(pars[1], pars[2], pars[3], pars[4]))
      
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
      
      # Save intermediate results after each iteration
      res_amo[[parsetnr]][[n_ix+1]] <- list(
        win_df = win_df,
        na_logL = na_logL,
        na_medL = na_medL
      )
      saveRDS(res_amo, file = "res_amo.rds")
    }
    
    # Reset restart index for next sample size
    restart_i <- 1
  }
  
  # Reset restart sample size for next parameter set
  restart_n <- nvec[1]
}

#-------------------------------------------------------------------------------
# When simulations are done, save the final results
saveRDS(res_amo, file = "res_amo.rds")
