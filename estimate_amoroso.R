#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### Load functions
require(AmoRosoDistrib)
require(dplyr)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

##----------------------------------------------------
# Define Amoroso negative log likelihood (LL) function
#-----------------------------------------------------

get_negLL_amo <- function(data, params) {
  
  # Amoroso parameters
  a <- params[1]
  l <- params[2]
  c <- params[3]
  mu <- params[4]
  
  # Negative log-likelihood of the data under that Amoroso
  negLL <- sum(-log(dgg4(data, a, l, c, mu)))
  
  return(negLL)
}

#----------------------------
# Define Amoroso BIC function
#----------------------------

get_BIC_amoroso <- function(data, params) {
  
  N <- length(data)
  P <- length(params)
  LL <- -get_negLL_amo(data, params)
  BIC <- LL + log(N) * (P + 1)
  
  return(BIC)
}

#--------------------------
# Create console info texts
#--------------------------

access_info <- "\nHOW TO ACCESS FULL OUTPUT:\n
  To access all of the output of the estimate_amoroso() function, you have to
  assign the results to an object like this:\n
  res <- estimate_amoroso(my_vector)\n
  Then you can access the following:\n
  - res$all_models: tibble that contains all 18 estimated Amorosos (9 methods
  x 2 parameter spaces (- and +)\n
  - res$min_BIC_models: tibble that contains, for each estimation method, the
  Amoroso fit with the lower BIC (i.e., either in -ve or +ve parameter space)\n
  - res$max_L_models: tibble that contains, for each estimation method, the
  Amoroso fit with the higher likelihood (i.e., either in -ve or +ve parameter
  space)\n
  - res$min_BIC_model: single-row tibble that contains the Amoroso fit with the
  lowest BIC overall\n
  - res$max_L_model: single-row tibble that contains the Amoroso fit with the
  highest likelihood overall\n\n"

plot1_info <- "\nABOUT THE PLOT:\n
  The plot shows a histogram of the data (grey bars),
  the nonparametric R Kernel density estimator fit (dark grey line), (optionally)
  the Amoroso fits from the initial parameter estimates (as described by Combes
  et al. (2022) and the Amoroso fits from the various estimation methods described
  by Combes et al. (2022). Each coloured line corresponds to the Amoroso fit of
  one method. Since each method yields two fits to the data - one in +ve
  and the other in -ve parameter space), the plot shows only the one that fits
  the data better according to the selected criterion (default 'ML' for maximum
  likelihood, alternatively 'BIC').\n\n"

plot2_info <- "\nABOUT THE PLOTS:\n
  The grid contains the Amoroso fits of all estimation methods that resulted
  in a valid fit (i.e., no mathematical errors). Therefore, each plot shows
  either the one or two valid fit(s) from one method. The dark-grey line is the
  fit of the nonparametric R Kernel density.\n\n"

plot3_info <- "\nABOUT THE PLOT:\n
  The plot shows the best Amoroso overall, as determined by either the maxmimum
  likelihood (criterion = 'ML') or the lowest BIC (criterion = 'BIC').\n\n"


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

estimate_amoroso <- function(vec = NULL,
                             n = 512,
                             criterion = "ML",
                             plot = 0,
                             legend = "topright",
                             breaks = 20,
                             print_results = FALSE,
                             include.init = FALSE) {
  
  # criterion = 'ML': select best model according to maximum likelihood
  # criterion = 'BIC': select best model according to minimum BIC
  # -> 'criterion' only relevant when plot = 1 or plot = 3
  
  # plot = 0: no plot
  # plot = 1: one plot with best Amoroso per method (+ve or -ve par space)
  # plot = 2: 3x3 grid with one method per plot
  # plot = 3: one plot with best Amoroso overall
  
  # breaks = n; n is the number of breaks in the histogram of the data
  
  # include.init = TRUE/FALSE: compatible with plot = 1, you can chose to either
  # add the 2 Amorosos from the initial parameter estimates (TRUE) to the plot
  # or not (FALSE)
  
  
  ##############################################################
  ## Make function that can handle failing estimation methods ##
  ##############################################################
  
  fit_safely <- function(method, init_value) {
    tryCatch({
      do.call(paste0("fit.", method), list(x, init_value))
    }, error = function(e) {
      message("Failed to fit.", method, ": ", e$message)
      return(NA)
    })
  }
  
  
  #####################
  ## Define data (x) ##
  #####################
  
  ## If data supplied in correct format:
  
  if (is.vector(vec) && is.numeric(vec)) {
    
    # Remove any NAs
    x <- na.omit(vec)
    
    # If NAs removed, print warning
    num_nas_removed <- length(vec) - length(na.omit(vec))
    if(num_nas_removed > 0) {
      message(c("WARNING: ", as.character(num_nas_removed),
                " NAs were removed from the data.\n"))
      cat("-----------------------------------------------------------------\n")
    }
  
  ## Else STOP
    
  } else {
    stop("'vec' must be a numeric vector. \n\n")
  }
  
  ## Nr of observations (after any NAs removed)
  length_x <- length(x)
  
  
  ######################################################################
  ## Estimate parameters for MLE and MDE methods in +ve and -ve space ##
  ######################################################################
  
  ## Find initializing parameters
  
  # -> For a > 0
  init.a.pos = init.theta(data = x, -20, 20, length = 1000, a.pos = TRUE)
  init.pos = init.a.pos[1:4]
  # -> For a < 0
  init.a.neg = init.theta(data = x, -20, 20, length = 1000, a.pos = FALSE)
  init.neg = init.a.neg[1:4]
  
  # List of methods and initial values
  methods <- c("mle","mkle","mjse","mhe","mwe","msqe","mhdfe","mwdfe","msqdfe")
  init_values <- list(init.pos = init.pos, init.neg = init.neg)
  
  
  ## Make list with all Amoroso fits
  
  # -> if method succeeds, assign list
  # -> if methods fails, assign NA
  
  fit_list <- list(
    # MLE
    mleresp = fit_safely("mle", init.pos), 
    mleresn = fit_safely("mle", init.neg),
    
    # Kullback-Leibler (PDF)
    mkleresp = fit_safely("mkle", init.pos), 
    mkleresn = fit_safely("mkle", init.neg),
    
    # Jensen-Shanon (PDF)
    mjseresp = fit_safely("mjse", init.pos), 
    mjseresn = fit_safely("mjse", init.neg),
    
    # Hellinger (PDF)
    mheresp = fit_safely("mhe", init.pos), 
    mheresn = fit_safely("mhe", init.neg),
    
    # Wasserstein (PDF)
    mWassderesp = fit_safely("mwe", init.pos), 
    mWassderesn = fit_safely("mwe", init.neg),
    
    # Squared (PDF)
    msqeresp = fit_safely("msqe", init.pos), 
    msqeresn = fit_safely("msqe", init.neg),
    
    # Hellinger (CDF)
    mhecresp = fit_safely("mhdfe", init.pos), 
    mhecresn = fit_safely("mhdfe", init.neg),
    
    # Wasserstein (CDF)
    mwecresp = fit_safely("mwdfe", init.pos), 
    mwecresn = fit_safely("mwdfe", init.neg),
    
    # Squared (CDF)
    msqdferesp = fit_safely("msqdfe", init.pos), 
    msqdferesn = fit_safely("msqdfe", init.neg)
  )
  
  
  ###################################################
  ## Make tibble with all estimated Amoroso models ##
  ###################################################

  all_models_tib <- tibble(
    
    # Method names
    method = rep(c("MLE", "Kullback-Leibler PDF", "Jensen-Shanon PDF",
                   "Hellinger PDF", "Wasserstein PDF", "Squared PDF",
                   "Hellinger CDF", "Wasserstein CDF", "Squared CDF"),
                 each = 2),
    
    # Method names ID
    method_ID = rep(c("MLE", "KL-PDF", "JS-PDF", "HELL-PDF", "WASS-PDF",
                      "SQU-PDF", "HELL-CDF", "WASS-CDF", "SQU-CDF"),
                    each = 2),
    
    # Parameter space
    space = rep(c("+", "-"), times = 9)
    
    ) %>% mutate(
      
      # Add parameters, negative LL and BIC of each method fit
      a = sapply(fit_list, function(fit)
        if (is.list(fit) && "par" %in% names(fit)) fit$par[1] else NA),
      l = sapply(fit_list, function(fit)
        if (is.list(fit) && "par" %in% names(fit)) fit$par[2] else NA),
      c = sapply(fit_list, function(fit)
        if (is.list(fit) && "par" %in% names(fit)) fit$par[3] else NA),
      mu = sapply(fit_list, function(fit)
        if (is.list(fit) && "par" %in% names(fit)) fit$par[4] else NA),
      negLL = sapply(fit_list, function(fit) {
        if (is.list(fit) && "par" %in% names(fit)) {
          get_negLL_amo(x, c(fit$par))
        } else { NA }}),
      BIC = sapply(fit_list, function(fit) {
        if (is.list(fit) && "par" %in% names(fit)) {
          get_BIC_amoroso(x, fit$par)
        } else { NA }})
    )
  
  ## Remove methods with NA or Inf parameters
  
  # Save IDs
  na_par_rows <- all_models_tib %>%
    filter(is.na(a) | is.na(l) | is.na(c) | is.na(mu)) %>%
    pull(method_ID)
  inf_par_rows <- all_models_tib %>%
    filter(!is.finite(a) | !is.finite(l) | !is.finite(c) | !is.finite(mu)) %>%
    pull(method_ID)
  
  # Remove rows
  all_models_tib <- all_models_tib %>%
    filter(!(is.na(a) | is.na(l) | is.na(c) | is.na(mu))) %>%
    filter(is.finite(a) & is.finite(l) & is.finite(c) & is.finite(mu))
  
  # Print warning messages
  if (length(na_par_rows) > 0) {
    warning(paste("The following methods had NA parameters and their rows were",
                  "removed): \n", paste(na_par_rows, collapse = ", ")))
  }
  if (length(inf_par_rows) > 0) {
    warning(paste("The following methods had Inf or -Inf parameters and their",
                  "rows were removed):\n", paste(na_par_rows, collapse = ", ")))
  }
  
  ## Remove methods with failed likelihood computation
  
  # Save IDs
  nan_ll_rows <- all_models_tib %>%
    filter(is.nan(negLL)) %>%
    pull(method_ID)
  
  # Remove rows
  all_models_tib <- all_models_tib %>%
    filter(!is.nan(negLL))
  
  # Print warning message
  if (length(nan_ll_rows) > 0) {
    warning(paste("The following methods' likelihood calculation failed (their",
    "rows were removed): \n", paste(nan_ll_rows, collapse = ", ")))
  }
  
  
  ###################################################
  ## Make dataframe of lowest BIC model PER METHOD ##
  ###################################################
  
  min_BIC_models_tib <- all_models_tib %>%
    group_by(method) %>%
    # Get the row with the minimum BIC per method
    filter(BIC == min(BIC)) %>%
    ungroup()
  
  
  ##################################################
  ## Make dataframe of higher LL model PER METHOD ##
  ##################################################
  
  max_L_models_tib <- all_models_tib %>%
    # Filter out rows where parameters are NA
    filter(!is.na(a) & !is.na(l) & !is.na(c) & !is.na(mu)) %>%
    group_by(method) %>%
    # Get the row with the minimum negLL (maximum likelihood)
    filter(negLL == min(negLL)) %>%
    ungroup()
  
  
  ###############################################
  ## Make dataframe of lowest BIC model OF ALL ##
  ###############################################
  
  min_BIC_model_tib <- min_BIC_models_tib %>%
    filter(BIC == min(BIC))
  
  
  #########################################
  ## Make dataframe of highest LL OF ALL ##
  #########################################
  
  max_L_model_tib <- max_L_models_tib %>%
    filter(negLL == min(negLL))
  
  
  ###########
  ## PLOTS ##
  ###########
  
  # Define x
  xx <- seq(min(density(x)$x), max(density(x)$x), length = n)
  
  # Fit R KDE
  rdens <- density(x)
  
  # Use R KDE to define a ymax that allows for some space at the top
  ybuffer <- 0.3*(range(rdens$y)[2]-range(rdens$y)[1])
  ymax <- max(rdens$y)+ybuffer
  
  #---------#
  # NO PLOT #
  #---------#
  
  if (plot == 0) {
    
    if (print_results == TRUE) {
      cat("no plot",
          "---------------------------------------------------------------\n\n")
    }
    
  #------------------------------------------------
  # PLOT = 1: One plot with best Amoroso per method
  #------------------------------------------------
  
    # -> the 'best' per method is in either +ve or -ve par space
    
  } else if (plot == 1) {
    
    # Select best model either according to BIC or ML
    if(criterion == "BIC") {
      win_models_tib <- min_BIC_models_tib
    } else if (criterion == "ML") {
      win_models_tib <- max_L_models_tib
    } else {
      stop("criterion must be either 'BIC' or 'ML'")
    }
    
    # Print plot info in console
    cat(plot1_info,
        "-------------------------------------------------------------------\n")
    
    # Plot settings
    par(cex.main = 1.4, cex.axis = 1, cex.lab = 1.2, bty = "n", font.lab = 2)
    
    # Create title
    title <- paste0("Amoroso fits to '", deparse(substitute(vec)), "'")

    # Create histogram
    hist(x, xlim = c(min(xx),max(xx)), ylim = c(0,ymax),
         probability = TRUE, main = title, breaks = breaks, axes = FALSE,
         col = "grey95", border = "grey85")
    axis(1)
    axis(2, las = 1)
    
    # Add R KDE
    lines(rdens$x, rdens$y, col = "grey70", lty = 1, lwd = 2)
    
    
    
    # Optional: include initial estimate Amorosos
    if (include.init == FALSE) {
      legend_names <- c("R density()")
      legend_colors <- c("grey70")
      lty_vec <- c(1, rep(1, nrow(win_models_tib)))
    } else {
      points(xx, dgg4(xx, init.pos[1], init.pos[2], init.pos[3], init.pos[4]),
             type = "l", lwd = 2, col = "black", lty = 1)
      points(xx, dgg4(xx, init.neg[1], init.neg[2], init.neg[3], init.neg[4]),
             type = "l", lwd = 2, col = "black", lty = 2)
      legend_names <- c("Initial Est.(+)", "Initial Est.(-)", "R density()")
      legend_colors <- c("black", "black", "grey70")
      lty_vec <- c(1, 2, 1, rep(1, nrow(win_models_tib)))
    }
    
    # Create colors for Amoroso fits
    model_colors <- rainbow(nrow(win_models_tib))
    
    # Loop to add each winning model
    for (i in 1:nrow(win_models_tib)) {
      # Create name for legend
      name <- paste0(win_models_tib$method[i], " (",
                     win_models_tib$space[i], ")")
      # Get model parameters
      pars <- c(win_models_tib$a[i], win_models_tib$l[i],
                win_models_tib$c[i], win_models_tib$mu[i])
      # Plot the winning model with a unique color
      points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]), 
             type = "l", lwd = 2, col = model_colors[i], lty = 1)
      # Add the model name to the legend
      legend_names <- c(legend_names, name)
      legend_colors <- c(legend_colors, model_colors[i])
    }
    
    # Add legend
    legend(legend, legend = legend_names,
           lwd = 2,
           col = legend_colors,
           lty = lty_vec,
           bty = "n",
           cex = 0.7)
    
  #-----------------------------------------------------
  # PLOT = 2: 3x3 grid of plots with one plot per method
  #-----------------------------------------------------
    
  } else if (plot == 2) {
    
    # Print plot info
    if (print_results == TRUE) {
      cat(plot2_info,
          "-----------------------------------------------------------------\n")
    }
    
    # Make 3x3 grid
    par(mfrow = c(3, 3), cex.axis = 0.8)
    
    # Loop through methods
    for (i in unique(all_models_tib$method)) {
      
      # Get all fits that exist for that method (1 or 2)
      meth_tib <- all_models_tib %>% filter(method == i)
      
      # If the method has two valid fits:
      if (nrow(meth_tib) == 2) {
        
        # Extract parameters of model in +ve space
        par_pos <- meth_tib %>%
          filter(space == "+") %>%
          select(a, l, c, mu) %>%
          unlist(use.names = FALSE)
        
        # Extract parameters of model in -ve space
        par_neg <- meth_tib %>%
          filter(space == "-") %>%
          select(a, l, c, mu) %>%
          unlist(use.names = FALSE)
        
        # Make the method name the title
        title <- meth_tib$method[[1]]
        
        # Make histogram
        hist(x, xlim = c(min(xx), max(xx)), ylim = c(0,ymax),
             probability = TRUE, main = title, breaks = breaks, axes = FALSE,
             col = "grey95", border = "grey85")
        axis(1)
        axis(2, las = 1)
        
        # Add R KDE
        lines(rdens$x, rdens$y, col = "grey70", lty = 1, lwd = 1.5)
        
        # Add Amoroso in +ve space
        points(xx, dgg4(xx, par_pos[1], par_pos[2], par_pos[3], par_pos[4]),
               type = "l", lwd = 1.5, col = "sienna3", lty = 1)
        
        # Add Amoroso in -ve space
        points(xx, dgg4(xx, par_neg[1], par_neg[2], par_neg[3], par_neg[4]),
               type = "l", lwd = 1.5, col = "steelblue3", lty = 1)
        
        # Add legend
        legend(legend, legend = c("R density()", "α > 0", "α < 0"),
               col = c("grey70", "sienna3", "steelblue3"), lty = c(1, 1, 1),
               lwd = 1.5, bty = "n", cex = 0.8)
        
      
      # If the method has only valid fit
      } else {
        
        # Extract parameters of single Amoroso model
        par <- meth_tib %>%
          slice(1) %>%
          select(a, l, c, mu) %>%
          unlist(use.names = FALSE)
        
        # Make the method name the title
        title <- meth_tib$method[[1]]
        
        # Make histogram
        hist(x, xlim = c(min(xx), max(xx)), ylim = c(0,ymax),
             probability = TRUE, main = title, breaks = breaks, axes = FALSE,
             col = "grey95", border = "grey85")
        axis(1)
        axis(2, las = 1)
        
        # Add R KDE
        lines(rdens$x, rdens$y, col = "grey70", lty = 1, lwd = 1.5)
        
        # Select colour for Amoroso model (+ve or -ve)
        ifelse(meth_tib$space == "+", col <- "sienna3", col <- "steelblue3")
        
        # Add Amoroso fit
        points(xx, dgg4(xx, pos_par[1], pos_par[2], pos_par[3], pos_par[4]),
               type = "l", lwd = 1.5, col = col, lty = 1)
        
        # Make legend
        # -> for Amoroso in +ve space
        if (meth_tib$space == "+") {
          legend_labels <- c("R density()", "α > 0")
          legend_colors <- c("grey70", "sienna3")
        } else {
        # -> for Amoroso in -ve space
          legend_labels <- c("R density()", "α < 0")
          legend_colors <- c("grey70", "steelblue3")
        }
        
        # Add legend
        legend(legend, legend = legend_labels,
               col = legend_colors, lty = c(1, 1, 1),
               lwd = 1.5, bty = "n", cex = 0.8)
      }
    }
    
  #------------------------------------
  # PLOT = 3: Plot only winning Amoroso
  #------------------------------------
    
  } else if (plot == 3) {
    
    # Use BIC or ML to select best Amoroso
    if(criterion == "BIC") {
      win_model_tib <- min_BIC_model_tib
    } else if (criterion == "ML") {
      win_model_tib <- max_L_model_tib
    } else {
      stop("criterion must be either 'BIC' or 'ML'")
    }
    
    # Print plot info
    if (print_results == TRUE) {
      cat(plot3_info,
          "---------------------------------------------------------------\n\n")
    }
    
    # Set plotting parameters
    par(cex.main = 1.4, cex.axis = 1, cex.lab = 1.2, bty = "n", font.lab = 2)
    
    # Extract name of best Amoroso
    method_name <- paste0(win_model_tib[["method"]]," (",
                          win_model_tib[["space"]], ")")
    
    # Create title
      title <- paste(method_name,"Amoroso")
    
    # Make histogram
    hist(x, xlim = c(min(xx),max(xx)), ylim = c(0,ymax),
         probability = TRUE, main = title, breaks = breaks, axes = FALSE,
         col = "grey95", border = "grey85", bty = "n")
    axis(1)
    axis(2, las = 1)
    
    # Add R KDE
    lines(rdens$x, rdens$y, col = "grey70", lty = 1, lwd = 2)
    
    # Extract parameters of best Amoroso
    pars <- win_model_tib %>%
      slice(1) %>%
      select(a, l, c, mu) %>% 
      unlist(use.names = FALSE)
    
    # Add Amoroso
    points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]), 
           type = "l", lwd = 3, col = "mediumpurple3", lty = 1)
    
    # Add legend
    legend_names <- c("R density()", method_name)
    legend_colors <- c("grey70","mediumpurple3")
    legend(legend, legend = legend_names,
           lwd = c(2,4),
           col = legend_colors,
           lty = c(1, 1),
           bty = "n",
           cex = 1.2)
    
    #----------------------------------------------------
    # ELSE: Error message: inadmissible plotting argument
    #----------------------------------------------------
    
  } else {
    
    stop("plot must be one of 0, (no plot),1 (all methods in one plot), 2
    (one plot per method) or 3 (one plot with lowest BIC Amoroso")
    
  }
  
  
  ######################
  ## Print and Return ##
  ######################

    # Print info on how to access full output
    if(print_results == TRUE) {
      cat(access_info, "\n")  
    }
  
  # Return all win models and final best model
  return(invisible(list(
    all_models = all_models_tib,
    min_BIC_models = min_BIC_models_tib,
    max_L_models = max_L_models_tib,
    min_BIC_model = min_BIC_model_tib,
    max_L_model = max_L_model_tib, 
    x = xx
  )))
  
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## Test

# set.seed(55)
# data <- rnorm(50, 2, 0.7)
# 

#dat <- palmerpenguins::penguins$flipper_length_mm
#estimate_amoroso(dat, plot = 0, criterion = "ML")
# res$all_models

