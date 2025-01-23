###############################################################################
# Packages and helper functions
###############################################################################

# Source functions
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/predict_mnorm.R"))

# Load packages
require(caret)
require(essHist)

# Helper function that generates Amoroso predictions
generate_amo_predictions <- function(test, method_name, res, fold) {
  if (length(res$modlist[[method_name]]) > 1) {
    pars <- res$modlist[[method_name]]$pars
    pred_y <- dgg4(test, pars[1], pars[2], pars[3], pars[4])
    if (anyNA(pred_y)) {
      nr_na <- sum(is.na(pred_y))
      pred_y[is.na(pred_y)] <- 0
      warning(paste("Fold", fold, ":", nr_na, "NA in the", method_name,
                    "predictions were replaced with zero."))
    }
    return(pred_y)
  } else {
    return(NA)
  }
}

###############################################################################
# get_pp() f
###############################################################################

get_pp <- function(
    dat, # data vector
    k = 5, # nr of folds
    generating_amoroso = NULL, 
    generating_normal = NULL,
    seed = 125 # seed to make fold creation reproducible
) {
  
  dat <- rnorm(100)
  k <- 2
  generating_amoroso <- NULL
  generating_amoroso <- NULL
  seed <- 24
  
  cat("--------------------------------------------------\n")
  cat("k-fold cross-validation with k =", k, "\n")
  cat("--------------------------------------------------\n")
  
  #---------------------------------------------------------------------------
  # Check if generating distribution is specified correctly
  #---------------------------------------------------------------------------
  add_prop_amo <- FALSE
  add_prop_norm <- FALSE
  
  if (!is.null(generating_amoroso) && !is.null(generating_normal)) {
    stop("Supply EITHER Amoroso OR Normal as data-generating distribution.")
  }
  
  if (!is.null(generating_amoroso)) {
    if (length(generating_amoroso) == 4 && is.numeric(generating_amoroso)) {
      add_prop_amo <- TRUE
      genpar <- generating_amoroso
    } else stop("Amoroso parameters must be 4 numbers (a, l, c, mu)")
  } else if (!is.null(generating_normal)) {
    if (length(generating_normal) == 2 && is.numeric(generating_normal)) {
      add_prop_norm <- TRUE
      genpar <- generating_normal
    } else stop("Normal parameters must be 2 numbers (mean, sd)")
  }
  
  #------------------------
  # Split data into k folds
  #------------------------
  # Get min and max value in data
  xmin <- min(dat)
  xmax <- max(dat)
  # Remove them from the data
  dat_stripped <- dat[dat != xmin & dat != xmax]
  # -> we do this to ensure that the density estimates are guaranteed to cover
  # the range of the and no zero likelihoods occur (we can't take log(0))
  # Divide data without min and max into folds
  set.seed(seed)
  folds <- createFolds(dat_stripped, k = k, list = TRUE, returnTrain = TRUE)
  
  #------------------------------------------------
  # Initialize tibble to store results of all folds
  #------------------------------------------------
  likelihood_tib <- tibble(
    method = c("rdens","scKDE_2infplus","mnorm",
               "amo_mle", "amo_hell_cdf", "amo_hell_pdf")
  )
  
  #---------------------------------------------------------------
  # Initialize vectors to store models with 0 or NA in predictions
  #---------------------------------------------------------------
  zero_pred_models <- c()
  na_pred_models <- c()
  
  #-------------------
  # Loop through folds
  #-------------------
  for (fold in 1:k) {
    
    #------------------------------------------
    # Make train and test data for current fold
    #------------------------------------------
    cat("\nfold =", fold, "\n")
    # Extract train data for current fold (from k-1 folds)
    train_indices <- folds[[fold]]
    train <- dat_stripped[train_indices]
    # Add the min and max data value to the train data
    train <- c(xmin,train,xmax)
    # Extract test data (from the remaining fold)
    test <- dat_stripped[-train_indices]
    cat("n =", length(train), "in train set;",
        "n =", length(test), "in test set\n")
    
    #-----------------------------------
    # Fit P and NP methods to train data
    #-----------------------------------
    res <- estimate_methods(dat=train, hist = TRUE, amorosocrit = "ML")
    
    #-----------------------------------------------
    # Make a list to store predictions for test fold
    #-----------------------------------------------
    pred_list <- list(rdens = NULL,
                      scKDE_2infplus = NULL,
                      mnorm = NULL,
                      amo_mle = NULL,
                      amo_hell_cdf = NULL,
                      amo_hell_pdf = NULL)
    
    #--------------------------------
    # Generate parametric predictions
    #--------------------------------
    # Generate predictions
    pred_list$amo_mle <- generate_amo_predictions(
      test, "amo_mle", res, fold)
    pred_list$amo_hell_cdf <- generate_amo_predictions(
      test, "amo_hell_cdf", res, fold)
    pred_list$amo_hell_pdf <- generate_amo_predictions(
      test, "amo_hell_pdf", res, fold)
    
    # Generate predictions from mixed normal
    mn_mod <- res$modlist_valid$mnorm
    pred_list$mnorm <- predict_mnorm(test, mn_mod, plot=F)
    
    #-----------------------------------
    # Generate nonparametric predictions
    #-----------------------------------
    res_interp <- res$modlist_valid_interp
    npfun_list <- list(rdens = NULL, scKDE_2infplus = NULL)
    
    for (i in 1:length(npfun_list)) {
      if (length(res_interp[[i]]) > 1) {
        npfun_list[[i]] <- splinefun(res_interp[[i]]$x,
                                     res_interp[[i]]$y,
                                     method = "monoH.FC")
      } else {
        npfun_list[[i]] <- NA
      }
    }
    
    for (i in 1:(length(pred_list)-4)) {
      if(!is.na(npfun_list[i])) {
        predvec <- npfun_list[[i]](test)
        predvec[is.na(predvec)] <- 0
        pred_list[[i]] <- predvec
      } else {
        pred_list[[i]] <- NA
      }
    }
    
    #----------------------------------------------
    # Calculate likelihood measures of current fold
    #----------------------------------------------
    # If we know the true distribution is Amoroso or normal, calculate the
    # true y for the test set
    if (add_prop_amo) {
      true_pred <- dgg4(test, genpar[1], genpar[2], genpar[3], genpar[4])
      true_logL <- sum(log(true_pred))
      true_medL <- median(true_pred)
    } else if (add_prop_norm) {
      true_pred <- dnorm(test, genpar[1], genpar[2])
      true_logL <- sum(log(true_pred))
      true_medL <- median(true_pred)
    }
    
    # Initialize vector for likelihood measures of all methods in this fold
    logL_vec_fold <- numeric(length(pred_list))
    medL_vec_fold <- numeric(length(pred_list))
    
    # For current fold, calculate the two likelihood measures for each method
    for (i in seq_along(pred_list)) {
      pred <- pred_list[[i]]
      method_name <- names(pred_list)[i]
      
      # If predictions include NA, calculate neither likelilihood measure
      if(anyNA(pred)) {
        na_pred_models <- append(na_pred_models, method_name)
        warning(paste("Fold", fold, ": Skipping calculations for method:",
                      method_name, "; Predictions include NA."))
        logL_vec_fold[i] <- NA
        medL_vec_fold[i] <- NA
        # If predictions include 0, calculate only median likelihood
      } else if (any(pred == 0)) {
        zero_pred_models <- append(zero_pred_models, method_name)
        warning(paste("Fold", fold, ": Skipping log-likelihood for method:",
                      method_name, "; Predictions include 0."))
        logL_vec_fold[i] <- NA
        medL_vec_fold[i] <- median(pred)
        # Otherwise calculate both
      } else {
        logL_vec_fold[i] <- sum(log(pred))
        medL_vec_fold[i] <- median(pred)
      }
    }
    
    #--------------------------------------------------
    # Add likelihood measures of current fold to tibble
    #--------------------------------------------------
    fold_colname_logL <- paste0("logL_f", fold)
    fold_colname_logL_prop <- paste0("logL_prop_f", fold)
    fold_colname_medL <- paste0("medL_f", fold)
    fold_colname_medL_prop <- paste0("medL_prop_f", fold)
    
    likelihood_tib[[fold_colname_logL]] <- logL_vec_fold
    likelihood_tib[[fold_colname_medL]] <- medL_vec_fold
    
    if(xor(add_prop_norm,add_prop_amo)) {
      likelihood_tib[[fold_colname_logL_prop]] <- logL_vec_fold/true_logL
      likelihood_tib[[fold_colname_medL_prop]] <- medL_vec_fold/true_medL
    }
  }
  
  #--------------------------------------------------------
  # Across all folds, calculate average likelihood measures
  #--------------------------------------------------------
  if(xor(add_prop_norm,add_prop_amo)) {
    likelihood_tib <- likelihood_tib %>%
      mutate(
        logL_avg = rowMeans(across(starts_with("logL_f"))),
        logL_prop_avg = rowMeans(across(starts_with("logL_prop_f"))),
        medL_avg = rowMeans(across(starts_with("medL_f"))),
        medL_prop_avg = rowMeans(across(starts_with("medL_prop_f")))
      )
  } else {
    likelihood_tib <- likelihood_tib %>%
      mutate(
        logL_avg = rowMeans(across(starts_with("logL_f"))),
        medL_avg = rowMeans(across(starts_with("medL_f")))
      )
  }
  
  return(list(
    likelihood_tib = likelihood_tib,
    all_na_models = unique(na_pred_models),
    all_zero_models = unique(zero_pred_models)
  ))
}


#------
dat <- rnorm(100,175,7)
res <- get_pp(dat)

