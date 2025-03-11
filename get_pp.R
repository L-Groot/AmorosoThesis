source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
get_pp <- function(
    dat, method = "k-fold", k = 5, prop_train = 0.8,
    generating_amoroso = NULL, generating_normal = NULL, 
    generating_exgauss = NULL, generating_mnorm_2comp = NULL,
    seed = 125)
  {
  
  #-----------------------------------------------------------------------------
  # Parameters :
  # generating Normal: c(mean, sd)
  # generating Amoroso: c(a l, c, mu)
  # generating bimodal mixed Normal: c(weight1, mean1, sd1, mean2, sd2)
  # generating ex-Gaussian: nc(mu, sigma, nu)
  
  #-----------------------------------------------------------------------------
  # set.seed(134)
  # dat <- rnorm(100)
  #method = "loocv"
  # k = 5
  # prop_train = 0.8
  #generating_amoroso = NULL
  #generating_normal = NULL
  #generating_exgauss <- c(pars[1],pars[2],pars[3])
  #generating_mnorm_2comp <- NULL
  # seed = 125
  
  #-----------------------------------------------------------------------------
  validate_inputs <- function() {
    # Count the number of non-null generating distributions
    num_generating_dists <- sum(!is.null(c(
      generating_amoroso, generating_normal, 
      generating_exgauss, generating_mnorm_2comp
    )))
    
    # Ensure exactly one generating distribution is specified
    if (num_generating_dists != 1) {
      stop("Specify exactly ONE generating distribution.")
    }
    
    # Check if correct number of parameters are supplied
    if (!is.null(generating_amoroso) && !(length(generating_amoroso) == 4)) {
      stop("Amoroso parameters must be 4 numeric values (a, l, c, mu).")
    }
    if (!is.null(generating_normal) && !(length(generating_normal) == 2)) {
      stop("Normal parameters must be 2 numeric values (mean, sd).")
    }
    if (!is.null(generating_exgauss) && !(length(generating_exgauss) == 3)) {
      stop("Ex-Gaussian parameters must be 3 numeric values (mu, sigma, nu).")
    }
    if (!is.null(generating_mnorm_2comp) && !(length(generating_mnorm_2comp) == 5)) {
      stop("Mixed Normal parameters must be 5 numeric values (comp1_prop, mu_comp1, sd_comp1, mu_comp2, sd_comp2).")
    }
    
    # Assign parameters
    if (length(generating_amoroso)==4) {
      genpar <- generating_amoroso
      gendist <- "amoroso"
    } else if (length(generating_normal)==2) {
      genpar <- generating_normal
      gendist <- "normal"
    } else if (length(generating_exgauss)==3) {
      genpar <- generating_exgauss
      gendist <- "exgaussian"
    } else if (length(generating_mnorm_2comp)==5) {
      genpar <- generating_mnorm_2comp
      gendist <- "mnorm_2comp"
    } else {
      genpar <- NULL
      gendist <- NULL
    }
    
    # Return generating distribution type and parameters
    return(list(genpar = genpar, gendist = gendist))
  }
  
  #-----------------------------------------------------------------------------
  calculate_true_likelihood <- function(test, genpar, distribution) {
    
    if (distribution == "amoroso") {
      pred <- dgg4(test, genpar[1], genpar[2], genpar[3], genpar[4])
      
    } else if (distribution == "normal") {
      pred <- dnorm(test, genpar[1], genpar[2])
      
    } else if (distribution == "exgaussian") {
      pred <- dexGAUS(test, mu = genpar[1], sigma = genpar[2], nu = genpar[3])
      
    } else if (distribution == "mnorm_2comp") {
      pred1 <- dnorm(test, mean = genpar[2], sd = genpar[3]) * genpar[1]
      pred2 <- dnorm(test, mean = genpar[4], sd = genpar[5]) * (1-genpar[1])
      pred <- pred1 + pred2
      
    } else {
      return(list(logL = NA, medL = NA))
    }
    
    # Return true log likelihoods and median likelihoods
    return(list(logL = sum(log(pred)), medL = median(pred)))
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
    
    res <- estimate_methods(train)
    res_interp <- res$modlist_valid_interp
    
    pred_list <- list(
      mnorm = predict_mnorm(test, res$modlist$mnorm, plot = FALSE),
      #amo_mle = predict_amo(test, "amo_mle", res),
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
  get_pp_measures <- function(train, test, genpar, gendist, np_methods) {
    
    # Calculate predictions for each model
    pred_list <- fit_and_predict(train, test, np_methods)
    
    # Calculate likelihood of prediction(s) for each model
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
    pp_measures_df <- do.call(rbind, lapply(pred_likelihoods, as.data.frame))
    
    # If true dist is known:
    if (!is.null(gendist)) {
      # Calculate likelihood of the test data under data-generating distribution
      true_likelihoods <- calculate_true_likelihood(test, genpar, gendist)
      
      # Add proportion of "true" likelihood measures to dataframe
      pp_measures_df <- pp_measures_df %>%
        mutate(prop_tru_logL = logL / true_likelihoods$logL,
               prop_tru_medL = medL / true_likelihoods$medL)
      
      # Calculate MSE of predictions
      mse_vec <- sapply(pred_list, function(pred) {
        if (length(pred) != length(test)) {
          warning("Prediction and test data length mismatch; skipping MSE calculation.")
          return(NA)
        } else {
          return(mean((pred - test)^2, na.rm = TRUE))
        }
      })
      
      # Add MSE column to the dataframe
      pp_measures_df <- pp_measures_df %>%
        mutate(mse = mse_vec)
    }
    
    # Return pp_measures_df
    return(pp_measures_df)
  }
  
  #-----------------------------------------------------------------------------
  run_kfold <- function(dat, k = 5, genpar, gendist,
                        np_methods = c("rdens", "scKDE_2infplus")) {
    
    # Remove min from data
    dat_stripped <- dat[dat != min(dat)]
    #dat_stripped <- dat[dat != min(dat) & dat != max(dat)]
    
    # Set seed for reproducible fold splitting
    set.seed(seed)
    
    # Divide data without min into k folds
    folds <- createFolds(dat_stripped, k = k, list = TRUE, returnTrain = TRUE)
    
    # Initialize lists to store results
    logL_results <- list()
    medL_results <- list()
    logL_prop_results <- list()
    medL_prop_results <- list()
    mse_results <- list()
    
    # Perform k-fold CV
    for (fold in seq_len(k)) {
      cat("Processing fold", fold, "\n")
      
      # Make train set (add min to it!)
      train <- c(min(dat), dat_stripped[folds[[fold]]])
      #train <- c(min(dat), dat_stripped[folds[[fold]]], max(dat))
      # -> always add min and max to train set
      test <- dat_stripped[-folds[[fold]]]
      
      # Get predictions for the current fold's test data
      pp_measures_df <- get_pp_measures(train, test, genpar, gendist, np_methods)
      
      # Extract logL and medL for the fold
      logL_results[[fold]] <- pp_measures_df$logL
      medL_results[[fold]] <- pp_measures_df$medL
      # Extract logL and medL as proportion of true likelihood measures
      logL_prop_results[[fold]] <- pp_measures_df$prop_tru_logL
      medL_prop_results[[fold]] <- pp_measures_df$prop_tru_medL
      # Extract MSE
      mse_results[[fold]] <- pp_measures_df$mse
    }
    
    # Combine results into data frames
    logL_df <- as.data.frame(do.call(cbind, logL_results))
    medL_df <- as.data.frame(do.call(cbind, medL_results))

    # Rename columns to indicate folds
    colnames(logL_df) <- paste0("logL_f", seq_len(k))
    colnames(medL_df) <- paste0("medL_f", seq_len(k))
    
    # Add average columns for each measure
    logL_df$avg_logL <- rowMeans(logL_df)
    medL_df$avg_medL <- rowMeans(medL_df)
    
    # Combine both logL and medL into one data frame (if required)
    avg_pp_df <- data.frame(logL_avg = logL_df$avg_logL,
                                     med_avg = medL_df$avg_medL)
    rownames(avg_pp_df) <- rownames(pp_measures_df)
    
    # If true distribution is supplied
    if (!is.null(gendist)) {
      # -> Add proportion likelihood measures
      logL_prop_df <- as.data.frame(do.call(cbind, logL_prop_results))
      medL_prop_df <- as.data.frame(do.call(cbind, medL_prop_results))
      colnames(logL_prop_df) <- paste0("logL_prop_f", seq_len(k))
      colnames(medL_prop_df) <- paste0("medL_prop_f", seq_len(k))
      logL_prop_df$avg_logL_prop <- rowMeans(logL_prop_df)
      medL_prop_df$avg_medL_prop <- rowMeans(medL_prop_df)
      avg_pp_df$logL_prop_avg <- logL_prop_df$avg_logL_prop
      avg_pp_df$medL_prop_avg <- medL_prop_df$avg_medL_prop
      # -> Add MSE
      mse_df <- as.data.frame(do.call(cbind, mse_results))
      colnames(mse_df) <- paste0("mse_f", seq_len(k))
      mse_df$avg_mse <- rowMeans(mse_df)
      avg_pp_df$mse_avg <- mse_df$avg_mse
    }
    
    # Return dataframe with average likelihood measures across folds
    return(avg_pp_df)
  }
  
  #-----------------------------------------------------------------------------
  run_splithalf <- function(dat, prop_train, genpar, gendist,
                            np_methods = c("rdens", "scKDE_2infplus")) {
    
    set.seed(seed)
    
    # Create initial train/test split
    inTrain <- createDataPartition(dat, p = prop_train, list = FALSE)
    train <- dat[inTrain]
    test <- dat[-inTrain]
    
    # Ensure that min is in train
    dat_min <- min(dat)
    
    # If min not yet in train, move it from test to train
    if (!(dat_min %in% train)) {
      train <- c(train, dat_min)
      test <- test[test != dat_min]
    }
    
    # Fit methods to train set, make predictions for test set and calculate
    # likelihood measures of the predictions
    pp_df <- get_pp_measures(train, test, genpar, gendist, np_methods)
  
    # Return dataframe
    return(pp_df)
  }
  
  #-----------------------------------------------------------------------------
  run_loocv <- function(dat, genpar, gendist,
                        np_methods = c("rdens", "scKDE_2infplus")) {
    
    # Initialize lists to store results
    logL_results <- list()
    medL_results <- list()
    logL_prop_results <- list()
    medL_prop_results <- list()
    mse_results <- list()
    skipped_iterations <- c()
    
    # Loop through all n test observations LOOCV
    for (i in seq_along(dat)) {
      
      # Skip iterations where the test observation is the min
      if (dat[i] == min(dat)) {
        skipped_iterations <- c(skipped_iterations, i)
        next
      }
      
      cat("Processing iteration", i, "\n")
      
      # Make train and test sets
      test <- dat[i]
      train <- dat[-i]
      
      # Get predictions for the current test observation
      pp_df <- get_pp_measures(train, test, genpar, gendist, np_methods)
      
      # Extract predictive performance measures for current iterations models
      logL_results[[i]] <- pp_df$logL
      medL_results[[i]] <- pp_df$medL
      logL_prop_results[[i]] <- pp_df$prop_tru_logL
      medL_prop_results[[i]] <- pp_df$prop_tru_medL
      mse_results[[i]] <- pp_df$mse
    }
    
    # Combine results into data frames
    logL_df <- as.data.frame(do.call(cbind, logL_results))
    medL_df <- as.data.frame(do.call(cbind, medL_results))
    logL_prop_df <- as.data.frame(do.call(cbind, logL_prop_results))
    medL_prop_df <- as.data.frame(do.call(cbind, medL_prop_results))
    mse_df <- as.data.frame(do.call(cbind, mse_results))
    
    # Add average columns for each measure
    logL_df$avg_logL <- rowMeans(logL_df)
    medL_df$avg_medL <- rowMeans(medL_df)
    logL_prop_df$avg_logL_prop <- rowMeans(logL_prop_df)
    medL_prop_df$avg_medL_prop <- rowMeans(medL_prop_df)
    mse_df$avg_mse <- rowMeans(mse_df)
    
    # Combine both logL and medL into one data frame (if required)
    avg_pp_df <- data.frame(logL_avg = logL_df$avg_logL,
                            medL_avg = medL_df$avg_medL,
                            mse_avg= mse_df$avg_mse)
    
    # If data-generating distribution is known:
    if (length(logL_prop_df) > 1 && length(medL_prop_df) > 1) {
      # Add average likelihood proportion measures
      avg_pp_df$logL_prop_avg <- logL_prop_df$avg_logL_prop
      avg_pp_df$medL_prop_avg <- medL_prop_df$avg_medL_prop
      # Add average MSE
      avg_pp_df$mse_avg <- mse_df$avg_mse
    }
    
    # Add column with method names
    avg_pp_df %>%
      mutate(method = rownames(pp_df)) %>%
      select(method, everything())

    # Make the rownames the methods
    rownames(avg_pp_df) <- rownames(pp_df)
    
    # Return dataframe with average likelihood measures
    return(avg_pp_df)
  }
  
  #-----------------------------------------------------------------------------
  genpar <- validate_inputs()$genpar
  gendist <- validate_inputs()$gendist

  if (method == "k-fold") {
    res <- run_kfold(dat, k, genpar, gendist)
  } else if (method == "split-half") {
    res <- run_splithalf(dat, prop_train, genpar, gendist)
  } else if (method == "loocv") {
    res <- run_loocv(dat, genpar, gendist)
    #max_logL_method <-
    #min_mse_method <-
  } else {
    stop("Invalid method. Use 'k-fold','split-half' or 'loocv'.")
  }
  
  return(res)
  
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#set.seed(70)
#dat <- rnorm(30)
# dat <- rgg4(20, a=4,l=1,c=7,mu=0)
# 
# 
#res <- get_pp(dat, method = "k-fold", k=5)
# res <- get_pp(dat, method = "k-fold", k=5, generating_amoroso = c(4,1,7,0))
# res <- get_pp(dat, method = "k-fold", k=5, generating_normal = c(0,1))
# 
# res <- get_pp(dat, method = "loocv")
# res <- get_pp(dat, method = "loocv", generating_amoroso = c(4,1,7,0))
#res <- get_pp(dat, method = "loocv", generating_normal = c(0,1))
# 
# res <- get_pp(dat, method = "split-half", prop_train = 0.5)
#res <- get_pp(dat, method = "split-half", prop_train = 0.5,
#            generating_amoroso = c(4,1,7,0))
# res <- get_pp(dat, method = "split-half", prop_train = 0.5,
#               generating_normal = c(0,1))
