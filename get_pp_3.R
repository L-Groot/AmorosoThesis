source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/predict_mnorm.R"))

library(caret)

#-------------------------------------------------------------------------------
get_pp <- function(
    dat, method = "k-fold", k = 5, prop_train = 0.8,
    generating_amoroso = NULL, generating_normal = NULL, seed = 125
) {
  
  dat <- rnorm(100,6,2)
  method = "k-fold"
  k = 5
  prop_train = 0.5
  generating_amoroso = NULL
  generating_normal = c(6,2)
  seed = 125
  
  #-----------------------------------------------------------------------------
  validate_inputs <- function() {
    if (!is.null(generating_amoroso) && !is.null(generating_normal)) {
      stop("Specify EITHER Amoroso OR Normal as the generating distribution.")
    }
    if (!is.null(generating_amoroso) && !(length(generating_amoroso) == 4)) {
      stop("Amoroso parameters must be 4 numeric values (a, l, c, mu).")
    }
    if (!is.null(generating_normal) && !(length(generating_normal) == 2)) {
      stop("Normal parameters must be 2 numeric values (mean, sd).")
    }
    if (length(generating_amoroso)==4) {
      genpar <- generating_amoroso
      gendist <- "amoroso"
    }
    else if (length(generating_normal)==2) {
      genpar <- generating_normal
      gendist <- "normal"
    } else {
      genpar <- NULL
      gendist <- NULL
    }
    return(list(genpar = genpar, gendist = gendist))
  }
  
  #-----------------------------------------------------------------------------
  genpar <- validate_inputs()$genpar
  gendist <- validate_inputs()$gendist
  
  #-----------------------------------------------------------------------------
  calculate_true_likelihood <- function(test, parameters, distribution) {
    if (distribution == "amoroso") {
      pred <- dgg4(test, genpar[1], genpar[2], genpar[3], genpar[4])
    } else if (distribution == "normal") {
      pred <- dnorm(test, genpar[1], genpar[2])
    } else {
      return(list(logL = NA, medL = NA))
    }
    list(logL = sum(log(pred)), medL = median(pred))
  }
  
  #-----------------------------------------------------------------------------
  predict_amo <- function(test, method_name, res) {
    if (method_name %in% names(res$modlist_valid)) {
      pars <- res$modlist[[method_name]]$pars
      pred <- dgg4(test, pars[1], pars[2], pars[3], pars[4])
      pred[is.na(pred)] <- 0
      return(pred)
    }
    return(rep(NA, length(test)))
  }
  
  #----------------------------------------------------------------------------
  fit_and_predict <- function(train, test, np_methods) {
    
    res <- estimate_methods(train, hist = TRUE)
    res_interp <- res$modlist_valid_interp
    
    cat("mnorm predictions: ", predict_mnorm(test, res$modlist_valid$mnorm, plot = FALSE), "\n\n")
    
    pred_list <- list(
      mnorm = predict_mnorm(test, res$modlist_valid$mnorm, plot = FALSE),
      amo_mle = predict_amo(test, "amo_mle", res),
      amo_hell_cdf = predict_amo(test, "amo_hell_cdf", res),
      amo_hell_pdf = predict_amo(test, "amo_hell_pdf", res)
    )
    
    npfun_list <- setNames(as.list(rep(NA,length(np_methods))), np_methods)
    
    for (model in np_methods) {
      np_model <- res_interp[[model]]
      npfun_list[[model]] <- list()
      # If fit is valid estimate np function
      if (length(np_model) > 1) {
        npfun_list[[model]] <- splinefun(res_interp[[model]]$x,
                                     res_interp[[model]]$y,
                                     method = "monoH.FC")
      } else {
        npfun_list[[model]] <- NA
      }
      # If fit is valid, generate np predictions
      if(!is.na(npfun_list[model])) {
        predvec <- npfun_list[[model]](test)
        predvec[is.na(predvec)] <- 0 # replace NA predictions with zero
        pred_list[[model]] <- predvec
        
      } else { # if the fit is invalid:
        pred_list[[model]] <- NA # -> add NA
      }
    }
    return(pred_list)
  }
  
  #-----------------------------------------------------------------------------
  get_likelihoods <- function(train, test, genpar, gendist, np_methods) {
    
    # Calculate predictions from all models
    pred_list <- fit_and_predict(train, test, np_methods)
    
    # Calculate likelihood of the test under all models
    pred_likelihoods <- lapply(pred_list, function(pred) {
      if (anyNA(pred)) {
        warning(paste("Fold [x]", ": Skipping log-likelihood and median
                        likelihood calculations for method:", "[method name]", ";
                        Predictions include NA."))
        return(list(logL = NA, medL = NA))
      } else if (any(pred == 0)){
        warning(paste("Fold", "[x]", ": Skipping log-likelihood calculation for
                        method:", "[method_name]", "; Predictions include 0. Still
                        calculating median likelihood."))
        return(list(logL = NA, medL = median(pred)))
      } else {
        return(list(logL = sum(log(pred)), medL = median(pred)))
      }
    })
    
    # Put likelihood measures of predictions in a df
    pred_likelihoods <- do.call(rbind, lapply(pred_likelihoods, as.data.frame))
    
    # If true dist is known: add proportion of "true" likelihood measures
    if (!is.null(gendist)) {
      # Calculate likelihood of the test data under data-generating distribution
      true_likelihoods <- calculate_true_likelihood(test, genpar, gendist)
      # -> log likelihood
      # -> median likelihood
      pred_likelihoods <- pred_likelihoods %>%
        mutate(prop_tru_logL = logL/true_likelihoods$logL,
               prop_tru_medL = medL/true_likelihoods$medL)
    }
    
    print(pred_likelihoods)
    
    # Return pred_likelihoods
    return(pred_likelihoods)
  }
  
  #-----------------------------------------------------------------------------
  run_kfold <- function(dat, k, genpar, gendist,
                        np_methods = c("rdens", "scKDE_2infplus")) {
    
    # Remove min and max from data
    dat_stripped <- dat[dat != min(dat) & dat != max(dat)]
    
    # Set seed for reproducible fold splitting
    set.seed(seed)
    
    # Divide data without min and max into folds
    folds <- createFolds(dat_stripped, k = k, list = TRUE, returnTrain = TRUE)
    
    # Initialize lists to store results
    logL_results <- list()
    medL_results <- list()
    logL_prop_results <- list()
    medL_prop_results <- list()
    
    # Perform k-fold CV and manually add min and max to each training set
    for (fold in seq_len(k)) {
      cat("Processing fold", fold, "\n")
      
      # Prepare train and test sets
      train <- c(min(dat), dat_stripped[folds[[fold]]], max(dat))
      test <- dat_stripped[-folds[[fold]]]
      
      # Get predictions for the current fold's test data
      likelihoods <- get_likelihoods(train, test, genpar, gendist, np_methods)
      
      # Extract logL and medL for the fold
      logL_results[[fold]] <- likelihoods$logL
      medL_results[[fold]] <- likelihoods$medL
      # Extract logL and medL as proportion of true likelihood measures
      logL_prop_results[[fold]] <- likelihoods$prop_tru_logL
      medL_prop_results[[fold]] <- likelihoods$prop_tru_medL
    }
    
    # Combine results into data frames
    logL_df <- as.data.frame(do.call(cbind, logL_results))
    medL_df <- as.data.frame(do.call(cbind, medL_results))
    logL_prop_df <- as.data.frame(do.call(cbind, logL_prop_results))
    medL_prop_df <- as.data.frame(do.call(cbind, medL_prop_results))
    
    # Rename columns to indicate folds
    colnames(logL_df) <- paste0("logL_f", seq_len(k))
    colnames(medL_df) <- paste0("medL_f", seq_len(k))
    colnames(logL_prop_df) <- paste0("logL_prop_f", seq_len(k))
    colnames(medL_prop_df) <- paste0("medL_prop_f", seq_len(k))
    
    # Add average columns for each measure
    logL_df$avg_logL <- rowMeans(logL_df)
    medL_df$avg_medL <- rowMeans(medL_df)
    logL_prop_df$avg_logL_prop <- rowMeans(logL_prop_df)
    medL_prop_df$avg_medL_prop <- rowMeans(medL_prop_df)
    
    # Combine both logL and medL into one data frame (if required)
    avg_likelihoods_df <- data.frame(logL_avg = logL_df$avg_logL,
                                     med_avg = medL_df$avg_medL,
                                     logL_prop_avg = logL_prop_df$avg_logL_prop,
                                     medL_prop_avg = medL_prop_df$avg_medL_prop)
    rownames(avg_likelihoods_df) <- rownames(likelihoods)
    
    return(avg_likelihoods_df)
  }
  
  #-----------------------------------------------------------------------------
  run_splithalf <- function(dat, prop_train, genpar, gendist,
                            np_methods = c("rdens", "scKDE_2infplus")) {
    
    set.seed(seed)
    
    # Sample a train set until both min and max are in it
    while (TRUE) {
      inTrain <- createDataPartition(dat, p = prop_train, list = FALSE)
      train <- dat[inTrain]
      if (all(c(min(dat), max(dat)) %in% train)) break
    }
    
    # Make the remaining observations the test set
    test <- dat[-inTrain]
    
    # Fit methods to train set, make predictions for test set and calculate
    # likelihood measures of the predictions
    res_likelihoods <- get_likelihoods(train, test, genpar, gendist, np_methods)
    
    likelihoods_df <- data.frame(
      logL = res_likelihoods$logL,
      medL = res_likelihoods$medL,
      logL_prop = res_likelihoods$prop_tru_logL,
      medL_prop = res_likelihoods$prop_tru_medL
    )
    
    rownames(likelihoods_df) <- rownames(res_likelihoods)
    
    # Return dataframe
    return(likelihoods_df)
  }
  
  #-----------------------------------------------------------------------------
  run_loocv <- function(dat, genpar, gendist,
                        np_methods = c("rdens", "scKDE_2infplus")) {
    
    # Initialize lists to store results
    logL_results <- list()
    medL_results <- list()
    logL_prop_results <- list()
    medL_prop_results <- list()
    skipped_iterations <- c()
    
    # Loop through all n test observations LOOCV
    for (i in seq_along(dat)) {
      
      # Skip iterations where the test observation is either the min or max
      if (dat[i] == min(dat) || dat[i] == max(dat)) {
        skipped_iterations <- c(skipped_iterations, i)  # Note down skipped iteration
        next
      }
      
      cat("Processing iteration", i, "\n")
      
      # Make train and test sets
      test <- dat[i]
      cat("test: ", test)
      train <- dat[-i]
      
      # Get predictions for the current test observation
      likelihoods <- get_likelihoods(train, test, genpar, gendist, np_methods)
      
      # Extract logL and medL for the iteration
      logL_results[[i]] <- likelihoods$logL
      medL_results[[i]] <- likelihoods$medL
      # Extract logL and medL as proportion of true likelihood measures
      logL_prop_results[[i]] <- likelihoods$prop_tru_logL
      medL_prop_results[[i]] <- likelihoods$prop_tru_medL
    }
    
    # Combine results into data frames
    logL_df <- as.data.frame(do.call(cbind, logL_results))
    medL_df <- as.data.frame(do.call(cbind, medL_results))
    logL_prop_df <- as.data.frame(do.call(cbind, logL_prop_results))
    medL_prop_df <- as.data.frame(do.call(cbind, medL_prop_results))
    
    # Add average columns for each measure
    logL_df$avg_logL <- rowMeans(logL_df, na.rm = TRUE)
    medL_df$avg_medL <- rowMeans(medL_df, na.rm = TRUE)
    logL_prop_df$avg_logL_prop <- rowMeans(logL_prop_df, na.rm = TRUE)
    medL_prop_df$avg_medL_prop <- rowMeans(medL_prop_df, na.rm = TRUE)
    
    # Combine both logL and medL into one data frame (if required)
    avg_likelihoods_df <- data.frame(logL_avg = logL_df$avg_logL,
                                     med_avg = medL_df$avg_medL)
    
    # If proportions were calculated, add the two columns to df
    if (length(logL_prop_df) == 18 && length(medL_prop_df) == 18) {
      avg_likelihoods_df$logL_prop_avg <- logL_prop_df$avg_logL_prop
      avg_likelihoods_df$medL_prop_avg <- medL_prop_df$avg_medL_prop
    }
    
    # Make the rownames the methods
    rownames(avg_likelihoods_df) <- rownames(likelihoods)
    
    # Return dataframe with average likelihood measures
    return(avg_likelihoods_df)
  }
  
  #-----------------------------------------------------------------------------
  validate_inputs()
  genpar <- if (!is.null(generating_amoroso)) generating_amoroso else generating_normal
  gendist <- if (!is.null(generating_amoroso)) "amoroso" else "normal"
  
  #-----------------------------------------------------------------------------
  
  if (method == "k-fold") {
    run_kfold(dat, k, genpar, gendist)
  } else if (method == "split-half") {
    run_splithalf(dat, prop_train, genpar, gendist)
  } else {
    stop("Invalid method. Use 'k-fold' or 'split-half'.")
  }
}


#-----

