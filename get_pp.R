#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------- 
# Source functions
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/pred_mnorm.R"))

# Load packages
require(AmoRosoDistrib)
require(caret)
require(tidyverse)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

get_pp <- function(
    dat, # data vector
    method = "k-fold", # "k-fold" or "split-half"
    k = 5, # nr of folds (for method = "k-fold)
    prop_train = 0.8, # prop. data in train set (for method = "split-half")
    generating_amoroso = NULL, #if data-generating dist. is Amoroso, add
    # predictive performance as proportion of maximum possible
    seed = 125 # seed to make fold creation reproducible
) 
  
  # set.seed(93)
  # dat <- rnorm(50, mean = 30, sd = 7)
  # method = "k-fold"
  # generating_amoroso = NULL
  # k = 5
  # prop_train = 0.8
  # seed = 125
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

{
  
  # If we know a data-generating Amoroso, we can express the log likelihood of
  # test data as a proportion of the log likelihood of that test data under the
  # 'true' model
  
  if (!is.null(generating_amoroso)) {
    if (length(generating_amoroso) == 4 && is.numeric(generating_amoroso)) {
      add_prop <- TRUE
      genpar <- generating_amoroso
    } else {
      stop("Amoroso parameters must be 4 numbers (a,l,c,mu)\n")
    }
  } else {
    add_prop <- FALSE
  }
  
  
  #-----------------------------------------------------------------------------  
  #-----------------------------------------------------------------------------
  
  ###############
  ### K- FOLD ###
  ###############
  
  
  if (method == "k-fold") {
    
    cat("----------------------------------\n")
    cat("K-fold cross-validation with k =", k,"\n")
    cat("----------------------------------\n")
    
    #------------------------
    # Split data into k folds
    #------------------------
    # Get min and max value in data
    xmin <- min(dat)
    xmax <- max(dat)
    # Remove them from the data
    dat_stripped <-  dat[dat != xmin & dat != xmax]
    # Set seed for reproducible fold splitting
    set.seed(seed)
    # Divide data without min and max into folds
    folds <- createFolds(dat_stripped, k = k, list = TRUE, returnTrain = TRUE)
    
    #------------------------------
    # Initialize tibble for results
    #------------------------------
    likelihood_tib <- tibble(
      method = c("rdens","bern1", "bern2",
                 "scKDE_uni", "scKDE_2inf", "scKDE_2infplus",
                 "mnorm", "amo_mle", "amo_hell_cdf", "amo_hell_pdf")
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
      
      # Print fold number
      cat("\nfold =", fold, "\n")
      # Extract train data for current fold (from k-1 folds)
      train_indices <- folds[[fold]]
      train <- dat_stripped[train_indices]
      # Add the min and max data value to the train data
      # -> this prevents zero predictions
      train <- c(xmin,train,xmax)
      # Extract test data (from the remaining fold)
      test <- dat_stripped[-train_indices]
      # Print nr ob observations in train and test set
      cat("n =", length(train), "in train set;",
          "n =", length(test), "in test set\n")
      
      # Fit P and NP methods to train data
      res <- estimate_methods(dat=train, hist = TRUE, amorosocrit = "ML")
      
      # Get interpolated density values (same x range)
      xy_interpolated <- res$modlist_valid_interp
      
      # Helper function that generates predictions from 3 Amorosos
      generate_amo_predictions <- function(test, method_name, res, fold) {
        
        # If the Amoroso is valid:
        if (length(res$modlist[[method_name]]) > 1) {
          # Get parameters
          pars <- res$modlist[[method_name]]$pars
          # Generate predictions
          pred_y <- dgg4(test, pars[1], pars[2], pars[3], pars[4])
          # Replace NA predictions with 0
          if (anyNA(pred_y)) {
            nr_na <- sum(is.na(pred_y))
            pred_y[is.na(pred_y)] <- 0
            warning(paste("Fold", fold, ":", nr_na, "NA in the", method_name,
                          "predictions were replaced with zero."))
          }
          return(pred_y)
          
        # If Amoroso is not valid:  
        } else {
          return(NA) # Make predictions NA
        }
      }
      
      # Generate predictions from MLE, Hellinger CDF and Hellinger PDF Amoroso
      pred_y_amo_mle <- generate_amo_predictions(
        test, "amo_mle", res, fold)
      pred_y_amo_hell_cdf <- generate_amo_predictions(
        test, "amo_hell_cdf", res, fold)
      pred_y_amo_hell_pdf <- generate_amo_predictions(
        test, "amo_hell_pdf", res, fold)
      
      # Generate predictions from mixed normal
      mn_mod <- res$modlist_valid$mnorm # Extract MN model
      pred_mnorm <- pred_mnorm(test, mn_mod, plot=F)
      
      # Extract list wiht interpolated x and y values
      res_interp <- res$modlist_valid_interp
      
      # Make a list to store continuous functions for the non-parametric fits
      npfun_list <- list(rdens = NULL,
                         bern1 = NULL,
                         bern2 = NULL,
                         scKDE_uni = NULL,
                         scKDE_2inf = NULL,
                         scKDE_2infplus = NULL)
      
      # Use splinefun() to make a continuous function from the NP estimates
      for (i in 1:length(npfun_list)) { #the first 6 fits in res are the np fits
        if (length(res_interp[[i]]) > 1) { #if the fit is valid, estimate function
          npfun_list[[i]] <- splinefun(res_interp[[i]]$x,
                                       res_interp[[i]]$y,
                                       method = "monoH.FC")
        } else { #if the fit is not valid, put NA
          npfun_list[[i]] <- NA
        }
      }
      
      # Make a list to store predictions for test fold
      pred_list <- list(rdens = NULL,
                        bern1 = NULL,
                        bern2 = NULL,
                        scKDE_uni = NULL,
                        scKDE_2inf = NULL,
                        scKDE_2infplus = NULL,
                        mnorm = NULL,
                        amo_mle = NULL,
                        amo_hell_cdf = NULL,
                        amo_hell_pdf = NULL)
      
      # Add predictions from NP fits to pred_list
      for (i in 1:(length(pred_list)-4)) { # leave out parametric fits
        
        if(!is.na(npfun_list[i])) { # if the fit is valid:
          predvec <- npfun_list[[i]](test) # -> add predictions
          # replace NA predictions with zero
          predvec[is.na(predvec)] <- 0
          # add predictions of this method to pred_list
          pred_list[[i]] <- predvec
          
        } else { # if the fit is invalid:
          pred_list[[i]] <- NA # -> add NA
        }
      }
      
      # Add Amoroso and Mixed Normal predictions to pred_list
        pred_list$mnorm <- pred_mnorm
        pred_list$amo_mle <- pred_y_amo_mle
        pred_list$amo_hell_cdf <- pred_y_amo_hell_cdf
        pred_list$amo_hell_pdf <- pred_y_amo_hell_pdf
      
        
      #### Calculate likelihood of this test folds' observations ####
      
      # Calculate likelihood of test set under data-generating Amoroso
      if (add_prop == TRUE) {
        
        true_pred <- dgg4(test, genpar[1], genpar[2], genpar[3], genpar[4])
        max_logL <- sum(log(true_pred))
        max_medL <- median(true_pred)
        
      }
      
      # Empty vector to store log-likelihoods for current fold
      logL_vec_fold <- numeric(length(pred_list))
      # Empty vector to store med-likelihoods for current fold
      medL_vec_fold <- numeric(length(pred_list))
      
      # Calculate log- and med-likelihood for current fold
      # -> log-likelihood is sum of log-likelihoods of observations in test fold
      # -> med-likelihood is median of the likelihoods of obs. in test fold
      for (i in seq_along(pred_list)) {
        pred <- pred_list[[i]]
        method_name <- names(pred_list)[i]
        # Check for NAs
        if(anyNA(pred)) {
          na_pred_models <- append(na_pred_models, method_name)
          warning(paste("Fold", fold, ": Skipping log-likelihood and median
                        likelihood calculations for method:", method_name, ";
                        Predictions include NA."))
          # Assign NA to both likelihood measures
          logL_vec_fold[i] <- NA
          medL_vec_fold[i] <- NA
        } else if (any(pred == 0)) {
          zero_pred_models <- append(zero_pred_models, method_name)
          warning(paste("Fold", fold, ": Skipping log-likelihood calculation for
                        method:", method_name, "; Predictions include 0. Still
                        calculating median likelihood."))
          # Assign NA to log-L measure
          logL_vec_fold[i] <- NA
          # -> Median likelihood of all methods predictions for this fold
          medL_vec_fold[i] <- median(pred)
        } else {
          # -> Log likelihood of all methods predictions for this fold
          logL_vec_fold[i] <- sum(log(pred)) 
          # -> Median likelihood of all methods predictions for this fold
          medL_vec_fold[i] <- median(pred) 
        }
      }
      
      # Create column names for current fold
      fold_colname_logL <- paste0("logL_f", fold)
      fold_colname_logL_prop <- paste0("logL_prop_f", fold)
      fold_colname_medL <- paste0("medL_f", fold)
      fold_colname_medL_prop <- paste0("medL_prop_f", fold)
      
      # Add column of log-likelihoods of all methods for current fold
      likelihood_tib[[fold_colname_logL]] <- logL_vec_fold
      # Add column with median-likelihoods of all methods of current fold
      likelihood_tib[[fold_colname_medL]] <- medL_vec_fold
      # Transform to proportion of "true" likelihood (under data-gen. model)
      if(add_prop == TRUE) {
        likelihood_tib[[fold_colname_logL_prop]] <- logL_vec_fold/max_logL
        likelihood_tib[[fold_colname_medL_prop]] <- medL_vec_fold/max_medL
      }
      
    }
    
    #------------------------------------------------------
    # Calculate average of likelihood measures across folds
    #------------------------------------------------------
    if (add_prop == TRUE) {
      
      # -> Likelihood measures expressed as proportion of "true" likelihoods
      
      # Add columns with average likelihood measures (across folds)
      likelihood_tib_full <- likelihood_tib %>%
        mutate(
          logL_avg = rowMeans(across(starts_with("logL_f"))),
          logL_prop_avg = rowMeans(across(starts_with("logL_prop_f"))),
          medL_avg = rowMeans(across(starts_with("medL_f"))),
          medL_prop_avg = rowMeans(across(starts_with("medL_prop_f")))
        )
      
      # Make new tibble with only the average columns
      likelihood_tib_avg <- likelihood_tib_full %>%
        select(method,logL_avg, logL_prop_avg, medL_avg, medL_prop_avg)
      
    } else {
      
      # -> Raw likelihood measures only (no proportions)
      
      # Add columns with average likelihood measures (across folds)
      likelihood_tib_full <- likelihood_tib %>%
        mutate(logL_avg = rowMeans(across(starts_with("logL_f")))) %>%
        mutate(medL_avg = rowMeans(across(starts_with("medL_f"))))
      
      # Make new tibble with only the average columns
      likelihood_tib_avg <- likelihood_tib_full %>%
        select(method,logL_avg, medL_avg)
    }
    
    
    #---------------
    # Return tibbles
    #---------------
    return(list(
      na_pred_models = unique(na_pred_models),
      zero_pred_models = unique(zero_pred_models),
      likelihood_tib_full = likelihood_tib_full,
      likelihood_tib_avg = likelihood_tib_avg))
    
    
  } else {
    
    
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    
    ##################
    ### SPLIT-HALF ###
    ##################
    
    cat("--------------------------------------------------------------\n")
    cat("Split-half CV (Proportion train = :", prop_train,"\n")
    cat("--------------------------------------------------------------\n")
    
    #-----------------------------
    # Split data in train and test
    #-----------------------------
    # Get min and max value in data
    xmin <- min(dat)
    xmax <- max(dat)
    # Remove them from the data
    dat_stripped <-  dat[dat != xmin & dat != xmax]
    # Set seed for reproducible data splitting
    set.seed(seed)
    # Split the data without the min and max into train and test
    inTrain <- createDataPartition(
      y = dat_stripped,
      p = prop_train, # prop. of data in the training set
      list = FALSE
    )
    # Split data in train and test
    train <- c(dat[inTrain],xmin,xmax)
    test <- dat[-inTrain]
    # Then add back the min and max to the train set
    # -> prevents zero predictions
    train <- c(train, xmin, xmax)
    
    #--------------------------------------------------------------------
    # Calculate log likelihood of test data under data-generating Amoroso
    #--------------------------------------------------------------------
    if (add_prop == TRUE) {
      true_pred <- dgg4(test, genpar[1], genpar[2], genpar[3], genpar[4])
      max_logL <- sum(log(true_pred))
      max_medL <- median(true_pred)
    }
    
    #---------------------------------------------------------------
    # Initialize vectors to store models with 0 or NA in predictions
    #---------------------------------------------------------------
    zero_pred_models <- c()
    na_pred_models <- c()
    
    #--------------------------------------------
    # Fit Amoroso and NP methods on training data
    #--------------------------------------------
    res <- estimate_methods(train, hist = TRUE, breaks = 20)
    
    #---------------------
    # Generate predictions
    #---------------------
    # Helper function that generates predictions from 3 Amorosos
    generate_amo_predictions <- function(test, method_name, res) {
      
      if (length(res$modlist_valid[[method_name]]) > 1) {
        
        # Get parameters
        pars <- res$modlist_valid[[method_name]]$pars
        # Generate predictions
        pred_y <- dgg4(test, pars[1], pars[2], pars[3], pars[4])
        # Replace NA predictions with 0
        if (anyNA(pred_y)) {
          nr_na <- sum(is.na(pred_y))
          pred_y[is.na(pred_y)] <- 0
          warning(paste(nr_na, "NA in the", method_name,
                        "predictions were replaced with zero."))
        }
        return(pred_y)
        
      } else {
        
        return(NA) # Make predictions NA
        
      }
    }
    
    # Generate predictions from mixed normal
    mn_mod <- res$modlist_valid$mnorm # Extract MN model
    pred_mnorm <- pred_mnorm(test, mn_mod, plot=F)
    
    # Generate predictions from Amorosos (MLE, Hellinger CDF, Hellinger PDF)
    pred_y_amo_mle <- generate_amo_predictions(test, "amo_mle", res)
    pred_y_amo_hell_cdf <- generate_amo_predictions(test, "amo_hell_cdf", res)
    pred_y_amo_hell_pdf <- generate_amo_predictions(test, "amo_hell_pdf", res)
    
    # Get interpolated Amoroso density values
    res_interp <- res$modlist_valid_interp
    
    # Make a list to store continuous functions for the non-parametric fits
    npfun_list <- list(rdens = NULL,
                       bern1 = NULL,
                       bern2 = NULL,
                       scKDE_uni = NULL,
                       scKDE_2inf = NULL,
                       scKDE_2infplus = NULL)
    
    # Use splinefun() to make a continuous function from the NP estimates
    for (i in 1:length(npfun_list)) { #the first 6 fits in res are the np fits
      if (length(res_interp[[i]]) > 1) { #if the fit is valid, estimate function
        npfun_list[[i]] <- splinefun(res_interp[[i]]$x,
                                     res_interp[[i]]$y,
                                     method = "monoH.FC")
      } else { #if the fit is not valid, put NA
        npfun_list[[i]] <- NA
      }
    }
    
    # Make a list to store predictions for test set
    pred_list <- list(rdens = NULL,
                      bern1 = NULL,
                      bern2 = NULL,
                      scKDE_uni = NULL,
                      scKDE_2inf = NULL,
                      scKDE_2infplus = NULL,
                      mnorm = NULL,
                      amo_mle = NULL,
                      amo_hell_cdf = NULL,
                      amo_hell_pdf = NULL)
    
    # Make predicitions from NP methods and add to pred_list
    for (i in 1:(length(pred_list)-4)) { # leave out parametric methods
      
      if(!is.na(npfun_list[i])) { # if the fit is valid:
        predvec <- npfun_list[[i]](test) # -> add predictions
        # replace NA predictions with 0
        predvec[is.na(predvec)] <- 0
        pred_list[[i]] <- predvec
        
      } else { # if the fit is invalid:
        pred_list[[i]] <- NA # -> add NA
      }
    }
    
    # Add Amoroso and mixed normal predictions to pred_list
    pred_list$mnorm <- pred_mnorm
    pred_list$amo_mle <- pred_y_amo_mle
    pred_list$amo_hell_cdf <- pred_y_amo_hell_cdf
    pred_list$amo_hell_pdf <- pred_y_amo_hell_pdf

    
    # Remove models that have NA predictions
    pred_list <- pred_list[!sapply(pred_list, function(x) any(is.na(x)))]
    
    #----------------------------------------------------------
    # For each model, calculate likelihood of test observations
    #----------------------------------------------------------
    # Empty vector to store log-likelihoods of test set OF EACH VALID MODEL
    logL_vec_testset <- numeric(length(pred_list))
    # Empty vector to store med-likelihoods of test set OF EACH VALID MODEL
    medL_vec_testset <- numeric(length(pred_list))
    
    # Loop through prediction vectors of all methods
    for (i in seq_along(pred_list)) {
      
      pred <- pred_list[[i]]
      method_name <- names(pred_list)[i]
      
      # Calculate likelihood measures for each method
      # Check for NAs
      if (any(pred == 0)) {
        zero_pred_models <- c(zero_pred_models, method_name)
        warning(paste("Skipping log-likelihood calculation for method:",
                      method_name, "; Predictions include 0. Still calculating
                      median likelihood."))
        # Assign NA to both likelihood measures
        logL_vec_testset[i] <- NA
        medL_vec_testset[i] <- median(pred)
      } else if (anyNA(pred)) {
        na_pred_models <- append(na_pred_models)
        warning(paste("Skipping log-likelihood and med-likelihood calculations
                      for method:", method_name, "; Predictions include NA."))
        # Assign NA for log-L measure
        logL_vec_testset[i] <- NA
        # -> Median likelihood of all methods predictions for the test set
        medL_vec_testset[i] <- NA
      } else {
        # -> Log likelihood of all methods predictions for the test set
        logL_vec_testset[i] <- sum(log(pred))
        # -> Median likelihood of all methods predictions for the test set
        medL_vec_testset[i] <- median(pred) 
      }
    }
    
    # Make tibble with likelihood measures
    likelihood_tib <- tibble(
      method = names(pred_list),
      logL_testset = logL_vec_testset,
      medL_testset = medL_vec_testset
    )

    # Express log-likelihood as proportions
    if(add_prop == TRUE) {
      likelihood_tib[["logL_prop_testset"]] <- logL_vec_testset/max_logL
      likelihood_tib[["medL_prop_testset"]] <- medL_vec_testset/max_medL
    }
    
    #---------
    # Returns
    #---------
    invisible(list(
      na_pred_models = na_pred_models,
      zero_pred_models = zero_pred_models,
      likelihood_tib = likelihood_tib
    ))
  }
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### Example usage ###

mydata <- rgg4(100, a=4, l=1, c=7, mu=0)
#(mydata <- rgg4(20, a=-1, l=2, c=2, mu=10))

#res <- get_pp(mydata, k=5)
res <- get_pp(mydata, k=5, generating_amoroso = c(-4,1,7,0))
#res <- get_pp(mydata, method="split-half")
#res <- get_pp(mydata, method = "split-half", generating_amoroso = c(-1,2,2,10))

#res$likelihood_tib_full
#res$likelihood_tib_avg
#res$likelihood_tib
  
# parset <- c(-4,1,-0.9,0)
# n <- 25
# dat <- rgg4(n, parset[1], parset[2], parset[3], parset[4])
# get_pp(dat)
