################################################################################
### Simulations: Data-generating ex-Gaussian
################################################################################

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in the list of results with what he have already
res_exgauss <- readRDS("res_exgauss_2.rds")

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
#nvec <- c(25,50)
nvec <- c(25,50,100,200)

# Nr of iterations per sample size and distribution
#nrep <- 3
nrep <- 100

# Make vector of seeds
seedvec <- seq(1,nrep)
seedvec[29] <- 177
seedvec[77] <- 23

# Initialize empty list for results
# res_exgauss <- list(
#   pars1 = list(
#     pars = c(),
#     n25 = list(),
#     n50 = list(),
#     n100 = list(),
#     n200 = list()
#   ),
#   pars2 = list(
#     pars = c(),
#     n25 = list(),
#     n50 = list(),
#     n100 = list(),
#     n200 = list()
#   ),
#   pars3 = list(
#     pars = c(),
#     n25 = list(),
#     n50 = list(),
#     n100 = list(),
#     n200 = list()
#   )
# )


# Define dynamic start points
#start_parsetnr <- 1  # Choose the parameter set number (1, 2, or 3)
#start_n <- 200      # Choose the sample size (25, 50, 100, 200)
#start_i <- 83        # Choose the iteration number (1 to nrep)

#-------------------------------------------------------------------------------
# Loop through parameter sets (starting at the chosen parset number)
# for (parsetnr in start_parsetnr:length(pars_list)) {
# 
#   # Get parameters
#   pars <- pars_list[[parsetnr]]
# 
#   # Store parameters in results
#   res_exgauss[[parsetnr]]$pars <- pars
# 
#   # Loop through sample sizes (starting at the chosen sample size)
#   for (n_ix in which(nvec == start_n):length(nvec)) {
# 
#     # Get current sample size
#     n <- nvec[n_ix]
#     print(n)
# 
#     # Initialize dataframes for results if they don't exist yet
#     if (length(res_exgauss[[parsetnr]][[n_ix + 1]])==0) {
#       res_exgauss[[parsetnr]][[n_ix + 1]] <- list(
#         win_df = data.frame(i = 1:nrep, max_logL_method = character(nrep), min_MSE_method = character(nrep)),
#         na_medL = data.frame(i = 1:nrep, rdens = numeric(nrep), scKDE_2infplus = numeric(nrep),
#                              mnorm = numeric(nrep), amo_hell_cdf = numeric(nrep), amo_hell_pdf = numeric(nrep)),
#         na_logL = data.frame(i = 1:nrep, scKDE_2infplus = numeric(nrep), mnorm = numeric(nrep),
#                              amo_hell_cdf = numeric(nrep), amo_hell_pdf = numeric(nrep))
#       )
#     }
# 
#     # Extract result storage
#     win_df <- res_exgauss[[parsetnr]][[n_ix + 1]]$win_df
#     na_medL <- res_exgauss[[parsetnr]][[n_ix + 1]]$na_medL
#     na_logL <- res_exgauss[[parsetnr]][[n_ix + 1]]$na_logL
# 
#     # Loop through iterations (starting at the chosen iteration)
#     for (i in if (parsetnr == start_parsetnr && n == start_n) start_i:nrep else 1:nrep) {
# 
#       # Print info of current iteration
#       cat("mu = ", pars[1], ", sigma = ", pars[2], ", nu = ", pars[3], " | ",
#           "n = ", n, " | ", "i = ", i, "\n")
# 
#       # Set seed dynamically
#       if (parsetnr == 1 && n == 100 && i == 45) {
#         set.seed(66)
#       } else if (parsetnr == 1 && n == 200 && i == 83){
#         set.seed(55)
#       } else {
#         set.seed(seedvec[i])
#       }
# 
#       # Simulate data
#       dat <- rexGAUS(n, mu = pars[1], sigma = pars[2], nu = pars[3])
# 
#       # Estimate methods and get PP measures
#       res <- get_pp(dat, method = "loocv",
#                     generating_exgauss = c(pars[1], pars[2], pars[3]))
# 
#       # Identify max-L and min-MSE method
#       max_logl_meth <- rownames(res)[which.max(ifelse(is.na(res$logL_avg), -Inf, res$logL_avg))]
#       min_mse_meth <- rownames(res)[which.min(ifelse(is.na(res$mse_avg), Inf, res$mse_avg))]
# 
#       # Store results in the corresponding dataframe
#       win_df$max_logL_method[i] <- max_logl_meth
#       win_df$min_MSE_method[i] <- min_mse_meth
# 
#       # Check if any method has NA for likelihood measures
#       methods <- rownames(res)
#       na_logL_flags <- as.numeric(is.na(res$logL_avg))
#       na_medL_flags <- as.numeric(is.na(res$medL_avg))
# 
#       # Store the NA flags in results dataframe
#       method_match <- intersect(methods, colnames(na_logL))
#       na_logL[i, method_match] <- na_logL_flags
#       na_medL[i, method_match] <- na_medL_flags
# 
#       # Save results **after every iteration**
#       res_exgauss[[parsetnr]][[n_ix + 1]] <- list(
#         win_df = win_df,
#         na_logL = na_logL,
#         na_medL = na_medL
#       )
#       saveRDS(res_exgauss, file = "res_exgauss.rds")
#     }
#   }
# }

#-------------------------------------------------------------------------------
# When simulations are done, save the list!
#saveRDS(res_exgauss, file = "res_exgauss.rds")

#-------------------------------------------------------------------------------

parsetnr <- 3

# Get parameters
pars <- pars_list[[parsetnr]]

# Store parameters in results
res_exgauss[[parsetnr]]$pars <- pars

# Loop through sample sizes (starting at the chosen sample size)
for (n_ix in 3:4) {
  
  # Get current sample size
  n <- nvec[n_ix]
  print(n)
  
  # Initialize dataframes for results if they don't exist yet
  if (length(res_exgauss[[parsetnr]][[n_ix + 1]])==0) {
    res_exgauss[[parsetnr]][[n_ix + 1]] <- list(
      win_df = data.frame(i = 1:nrep, max_logL_method = character(nrep), min_MSE_method = character(nrep)),
      na_medL = data.frame(i = 1:nrep, rdens = numeric(nrep), scKDE_2infplus = numeric(nrep),
                           mnorm = numeric(nrep), amo_hell_cdf = numeric(nrep), amo_hell_pdf = numeric(nrep)),
      na_logL = data.frame(i = 1:nrep, scKDE_2infplus = numeric(nrep), mnorm = numeric(nrep),
                           amo_hell_cdf = numeric(nrep), amo_hell_pdf = numeric(nrep))
    )
  }
  
  # Extract result storage
  win_df <- res_exgauss[[parsetnr]][[n_ix + 1]]$win_df
  na_medL <- res_exgauss[[parsetnr]][[n_ix + 1]]$na_medL
  na_logL <- res_exgauss[[parsetnr]][[n_ix + 1]]$na_logL
  
  start_ix <- ifelse(n == 100, 26, 1)
  
  # Loop through iterations (starting at the chosen iteration)
  for (i in start_ix:nrep) {
    
    # Print info of current iteration
    cat("mu = ", pars[1], ", sigma = ", pars[2], ", nu = ", pars[3], " | ",
        "n = ", n, " | ", "i = ", i, "\n")
    
    # Set seed dynamically
    # if (parsetnr == 1 && n == 100 && i == 45) {
    #   set.seed(66)
    # } else if (parsetnr == 1 && n == 200 && i == 83){
    #   set.seed(55)
    # } else if (parsetnr == 2 && n == 50 && i == 95) {
    #   set.seed(44)
    # } else if (parsetnr == 2 && n == 100 && i == 66){
    #   set.seed(32)
    # } else if (parsetnr == 3 && n == 25 && i == 82){
    #   set.seed(177)
    # } else if (parsetnr == 3 && n == 50 && i == 67){
    #   set.seed(43)
    # 
    
    if (parsetnr == 3 && n == 100 && i == 26){
      set.seed(44)
    } else {
      set.seed(seedvec[i])
    }
    
    # Simulate data
    dat <- rexGAUS(n, mu = pars[1], sigma = pars[2], nu = pars[3])
    
    # Estimate methods and get PP measures
    res <- get_pp(dat, method = "loocv",
                  generating_exgauss = c(pars[1], pars[2], pars[3]))
    
    # Identify max-L and min-MSE method
    max_logl_meth <- rownames(res)[which.max(ifelse(is.na(res$logL_avg), -Inf, res$logL_avg))]
    min_mse_meth <- rownames(res)[which.min(ifelse(is.na(res$mse_avg), Inf, res$mse_avg))]
    
    # Store results in the corresponding dataframe
    win_df$max_logL_method[i] <- max_logl_meth
    win_df$min_MSE_method[i] <- min_mse_meth
    
    # Check if any method has NA for likelihood measures
    methods <- rownames(res)
    na_logL_flags <- as.numeric(is.na(res$logL_avg))
    na_medL_flags <- as.numeric(is.na(res$medL_avg))
    
    # Store the NA flags in results dataframe
    method_match <- intersect(methods, colnames(na_logL))
    na_logL[i, method_match] <- na_logL_flags
    na_medL[i, method_match] <- na_medL_flags
    
    # Save results **after every iteration**
    res_exgauss[[parsetnr]][[n_ix + 1]] <- list(
      win_df = win_df,
      na_logL = na_logL,
      na_medL = na_medL
    )
    saveRDS(res_exgauss, file = "res_exgauss_2.rds")
  }
}

saveRDS(res_exgauss, file = "res_exgauss.rds")
