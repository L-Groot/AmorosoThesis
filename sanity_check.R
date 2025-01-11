# SANITY CHECK:

# When simulating data from Amoroso, ML-Amoroso should - on average - produce
# the maximum likelihood fit

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Source functions
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))

# Load packages
require(tidyverse)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### List of parameter sets
parsets <- list(
  # -> taken from proposal appendix
  # -> but removed any distributions that are not unimodal (exponential)
  # because for these the adjusted KDE methods and the Amoroso PDF will fail
  # -> adj KDE will keep searching for a solution that obeys the constraints
  # but not find one
  # The Amoroso CDFs can return a fit without a peak
  # The Bernstein and R density work but return a bad fit
  
  #list(c(-4,1,-0.9,0))
  #list(c(-1,1,2,10), c(-1,2,2,10))
  
  # a grows pos
  list(c(1,2,2,0), c(2,2,2,0), c(3,2,2,0), c(4,2,2,0)),
  
  # a grows neg
  list(c(-1,2,2,10), c(-2,2,2,10), c(-3,2,2,10), c(-4,2,2,10)),
  
  # l in a-
  list(
    #c(-1,0.1,2,10), c(-1,0.5,2,10),
    c(-1,1,2,10), c(-1,2,2,10),
    c(-1,10,2,10), c(-1,30,2,10)),
  
  # l in a+
  list(
    #c(1,0.1,2,0), c(1,0.5,2,0),
    c(1,1,2,0), c(1,2,2,0),
    c(1,10,2,0), c(1,30,2,0)),
  
  # c grows neg in a-
  list(
    #c(-4,1,-0.1,0),
    c(-4,1,-0.9,0), c(-4,1,-1.1,0), c(-4,1,-5,0),
    c(-4,1,-7,0), c(-4,1,-10,0)),
  
  # c grows pos in a-
  list(
    #c(-4,1,0.1,0),
    c(-4,1,0.9,0), c(-4,1,1.1,0), c(-4,1,5,0),
    c(-4,1,7,0), c(-4,1,10,0)),
  
  # c grows neg in a+
  list(
    #c(4,1,-0.1,0),
    c(4,1,-0.9,0), c(4,1,-1.1,0), c(4,1,-5,0),
    c(4,1,-7,0), c(4,1,-10,0)),
  
  # c grows pos in a+
  list(
    c(4,1,1.1,0), c(4,1,5,0),
    c(4,1,7,0), c(4,1,10,0)),
  
  # mu in a-
  list(c(-1,2,1,0), c(-1,2,1,1), c(-1,2,1,2), c(-1,2,1,3),
       c(-1,2,1,4), c(-1,2,1,5)),
  
  # mu in a+
  list(c(1,2,1,0), c(1,2,1,1), c(1,2,1,2), c(1,2,1,3),
       c(1,2,1,4), c(1,2,1,5))
)

### Sample sizes
nvec <- c(50, 75, 100, 150)

### Loop through different sample sizes and Amoroso parameter sets
for (n in nvec) {
  
  # Initialize an empty tibble for each sample size
  assign(paste0("results_tib_", n), tibble())
  
  for (parset in unlist(parsets, recursive =  FALSE)) {
    
    cat("N =", n,"\n")
    cat("PARS: a=", parset[1], ", l=", parset[2], ", c=", parset[3], ", mu=", parset[4], "\n")
    
    # Generate Amoroso data
    dat <- rgg4(n, a = parset[1], l = parset[2], c = parset[3], mu = parset[4])
    #print(dat)
    
    # Fit methods and assess predictive performance
    res <- get_pp(dat, generating_amoroso = parset)
    
    # Extract method with highest logL and method with highest medL
    best_methods <- res$likelihood_tib_avg %>%
      
      # Remove rows where any value is NA
      filter(!is.na(logL_avg), !is.na(logL_prop_avg), !is.na(medL_avg)) %>%
      
      summarise(
        # Highest logL_avg and corresponding method
        logL_1_meth = method[which.max(logL_prop_avg)],
        logL_1_prop = logL_prop_avg[which.max(logL_prop_avg)],
        
        # Second-highest logL_avg and corresponding method
        logL_2_meth = method[order(-logL_prop_avg)[2]],
        logL_2_prop = logL_prop_avg[order(-logL_prop_avg)[2]],
        
        # Third-highest logL_avg and corresponding method
        logL_3_meth = method[order(-logL_prop_avg)[3]],
        logL_3_prop = logL_prop_avg[order(-logL_prop_avg)[3]],
        
        # Highest medL_avg and corresponding method
        medL_1_meth = method[which.max(medL_avg)],
        medL_1 = medL_avg[which.max(medL_avg)],
        
        # Second-highest medL_avg and corresponding method
        medL_2_meth = method[order(-medL_avg)[2]],
        medL_2 = medL_avg[order(-medL_avg)[2]],
        
        # Third-highest medL_avg and corresponding method
        medL_3_meth = method[order(-medL_avg)[3]],
        medL_3_prop = medL_avg[order(-medL_avg)[3]]
      )
    
    # Add column with names of models that had NA and 0 predictions
    best_methods$na_pred <- list(res$na_pred_models)
    best_methods$zero_pred <- list(res$zero_pred_models)
    
    # Add parameters of data-generating Amoroso
    best_methods <- best_methods %>%
      mutate(
        a = parset[1],
        l = parset[2],
        c = parset[3],
        mu = parset[4]
      )  %>%
      relocate(a, l, c, mu) # make parameters the first four columns
    
    # Add results to tibble
    best_methods <- bind_rows(get(paste0("results_tib_", n)), best_methods)
    
    # Add the results as a new row to the tibble for the current sample size
    assign(paste0("results_tib_", n), best_methods)
  }
}

# Optionally, print each results tibble to verify
for (n in nvec) {
  cat("Results for N = ", n, ":\n")
  print(get(paste0("results_tib_", n)))
}


# Save all tibbles in a single .RData file
save(results_tib_50, results_tib_75, results_tib_100, results_tib_150, file = "results_tibbles.RData")
# Load the saved tibbles into the environment
#load("results_tibbles.RData")