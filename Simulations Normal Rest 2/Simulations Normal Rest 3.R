################################################################################
### Simulations: Data-generating Normal
################################################################################

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in the list of results with what he have already
res_normal <- readRDS("res_normal.rds")

# Load functions from Github
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))


# Define parameters for the normal distributions
par_norm_wide <- c(100, 15)
par_norm_medium <- c(100, 10)
par_norm_narrow <- c(100, 5)

# Put parameters in list
pars_list <- list(par_norm_wide, par_norm_medium, par_norm_narrow)

# Set nr of repetitions (per parameter set and sample size
#nrep <- 2
nrep <- 100

# Make seeds
seedvec <- seq(1,nrep)

# Initialize list for normal results
# res_normal <- list(
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

# Get parameters
pars <- c(100,15)

parsetnr <- 1

n_ix <- 4

# Store parameters in results
res_normal[[parsetnr]]$pars <- pars

# Get current sample size
n <- 200

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
  cat("mean = ", pars[1], "sd = ", pars[2], " | ", "n = ", n, " | ",
      "i = ", i, "\n")
  
  # Simulate data
  set.seed(seedvec[i])
  dat <- rnorm(n, mean = pars[1], sd = pars[2])
  
  # Estimate methods and get PP measures
  res <- get_pp(dat, method = "loocv", generating_normal = c(pars[1], pars[2]))
  
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
  
  # Update the results in the main list
  res_normal[[parsetnr]][[n_ix + 1]] <- list(win_df = win_df,
                                             na_logL = na_logL,
                                             na_medL = na_medL)
  
  # Save after every iteration
  saveRDS(res_normal, file = "res_normal.rds")
  
}


#-------------------------------------------------------------------------------
  
  
