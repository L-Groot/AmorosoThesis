# estimate_methods 3
# -> removed Bernstein
# -> removed adj. KDE "unimodal and adj KDE "twoInflections+"
# -> removed Amoroso MLE
# -> force to take only a>0 Amorosos

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Source functions
# -> for estimating Amoroso
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/mnorm_functions.R"))

# Load packages
# -> for Amoroso density function
require(AmoRosoDistrib)
# -> for adjusted KDE
require(scdensity) 
# -> for mixed normal estimation
require(mclust)
# -> for mixed normal density function
require(LaplacesDemon)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

estimate_methods <- function(dat = NULL,
                             plot = TRUE, hist = TRUE, breaks = 20,
                             amoinaplus = FALSE, #only consider amorosos in a>0?
                             minimal = FALSE,
                             plot_common_x = TRUE,
                             rug = FALSE,
                             main = NULL,
                             generatingnormal = NULL, #supply mean,sd)
                             generatingamoroso = NULL, #supply (a,l,c,mu)
                             xticks = NULL,
                             yticks = NULL
) {
  
  # # For testing
  # dat = rnorm(70)
  # plot = TRUE; hist = TRUE; breaks = 20
  # amoinaplus = TRUE
  # minimal = FALSE
  # plot_common_x = TRUE
  # rug = FALSE
  # main = NULL
  # generatingnormal = NULL # supply (mean,sd)
  # xticks = NULL
  
  ########################
  ### HELPER FUNCTION  ###
  ########################
  
  # Function that can handle errors in estimation
  safe_execute <- function(expr, object_name, data_vector) {
    tryCatch(
      {
        result <- eval(bquote(.(expr)), envir = list(dat = data_vector))
        return(result)
      },
      error = function(e) {
        cat(paste("Error with fitting", object_name, ":", e$message, ";\n",
                  "Other methods were still fit.\n"))
        return(NA) # assigns NA if method failed
      }
    )
  }
  
  ############################
  ### 1. ESTIMATE DENSITY  ###
  ############################
  
  # Remove NA
  num_nas_to_remove <- length(dat) - length(na.omit(dat))
  dat <- as.vector(na.omit(dat))
  
  # Print how many NAs were removed
  if(num_nas_to_remove > 0) {
    message(c("WARNING: ", as.character(num_nas_to_remove),
              " NAs were removed from the data.\n"))
    cat("--------------------------------------------------------------------\n")
  }
  
  # Get n
  n <- length(dat)
  
  #### Amoroso ####
  amo <- safe_execute(quote(
    estimate_amoroso(dat, plot=0, criterion="ML")), "amo", dat)
  amo_x <- amo$x
  if (amoinaplus) {
    # For each method, extract model in a>0 parameter space
    amo <- amo$all_models #extract df with all models
    amo <- amo %>% filter(space == "+") #extract only models in a>0 par space
  } else {
    # For each method, extract model with the higher likeliihood (a+ or a-)
    amo <- amo$max_L_models
  }

  
  #### Adjusted KDE ####
  scKDE_2infplus <- safe_execute(quote(
    scdensity(dat, constraint = "twoInflections+")), "scKDE_2infplus", dat)
  
  ##### R density ####
  rdens <- safe_execute(quote(
    density(dat)), "rdens", dat)
  
  ##### Mixed Normal #####
  mnorm <- safe_execute(quote(
    densityMclust(dat, plot=F)), "mnorm", dat)
  mnorm$x <- rdens$x
  mnorm$y <- predict_mnorm(mnorm$x, mnorm, plot=F)
  
  
  ############################
  ### 2. EXTRACT AMOROSOS  ###
  ############################
  
  # Helper function to extract a specific Amoroso from the big Amoroso list
  extract_amoroso <- function(amo, method_id) {
    if (method_id %in% amo$method_ID) {
      model <- amo %>% filter(method_ID == method_id)
      pars <- model %>%
        slice(1) %>%
        select(a, l, c, mu) %>%
        unlist(use.names = FALSE)
      name <- paste0(model$method, " (", model$space, ")")
      name_id <- paste0(model$method_ID, " (", model$space, ")")
      density_values <- dgg4(amo_x, pars[1], pars[2], pars[3], pars[4])
      density_values[is.na(density_values)] <- 0
      return(list(x = amo_x, y = density_values, pars = pars, method = name,
                  method_short = name_id))
    } else {
      warning(paste(method_id, "Amoroso absent"))
      return(NULL)
    }
  }
  
  # Extract Hellinger CDF and Hellinger PDF Amoroso in a>0 par space
  if (length(amo) > 1) { # If at least one Amoroso method produced a valid fit:
    if ("HELL-CDF" %in% amo$method_ID) {
      amo_hell_cdf <- extract_amoroso(amo, "HELL-CDF")
    } else {
      amo_hell_cdf <- NULL
    }
    if ("HELL-PDF" %in% amo$method_ID) {
      amo_hell_pdf <- extract_amoroso(amo, "HELL-PDF")
    } else {
      amo_hell_pdf <- NULL
    }
  }
  
  
  ############################
  ### 3. MAKE MODEL LISTS  ###
  ############################
  
  # List of all models
  modlist <- list(rdens = rdens,
                  scKDE_2infplus = scKDE_2infplus,
                  mnorm = mnorm,
                  amo_hell_cdf = amo_hell_cdf,
                  amo_hell_pdf = amo_hell_pdf
  )
  
  # List of only valid models
  # -> (i.e., the ones that didn't fail to fit)
  # -> We need this list to define the x range for the plots
  modlist_valid <- modlist[sapply(modlist, function(mod) length(mod) > 1)]
  
  # Get xmin and xmax across all valid models
  xmin <- min(sapply(modlist_valid, function(mod) min(mod$x)))
  xmax <- max(sapply(modlist_valid, function(mod) max(mod$x)))
  
  # Define x range
  xvals <- seq(xmin, xmax, length.out = 512)
  
  # Get ymin and ymax across all valid models
  ymaxes <- sort(sapply(modlist_valid, function(mod) max(mod$y)),
                 decreasing = TRUE)
  
  # If any of the valid models has a spike in the density estimate, cut that
  # fit off -> just make the max histogram value the ymax
  # if (ymaxes[1] > (3*ymaxes[length(ymaxes)])) {
  #   hist_list <- hist(dat, breaks=breaks, prob=TRUE)
  #   ymax <- max(hist_list$density) #-> make ymax the max histogram density
  #   warning("At least one density estimate has a spike which was cut off in the plots")
  # } else {
  #   ymax <- ymaxes[1]
  # }
  hist_list <- hist(dat, breaks=breaks, prob=TRUE)
  ymax <- max(hist_list$density)
  buffer <- 0.15*ymax
  ymax <- ymax + buffer
  
  # Interpolate the valid models so that they all cover the same x range
  modlist_valid_interp <- lapply(names(modlist_valid), function(name) {
    mod <- modlist_valid[[name]]
    ## For valid models: interpolate
    if (length(mod) > 1) {
      y <- approx(mod$x, mod$y, xout = xvals, rule = 1)$y
      y[is.na(y)] <- 0 # Replace NA density values with 0
      list(x = xvals, y = y)
      ## For invalid models: assign NA to model predictions
    } else {
      NA
    }
  })
  # Add model names
  names(modlist_valid_interp) <- names(modlist_valid)
  
  # Make xlim for plot
  xmin_dat <- min(dat)
  xmax_dat <- max(dat)
  xrange <- xmax_dat-xmin_dat
  xmin_plot <- xmin_dat-0.2*xrange
  xmax_plot <- xmax_dat+0.2*xrange
  
  #####################
  ### 4. MAKE PLOTS ###
  #####################
  
  if (plot == TRUE) {
    
    # Plot either common x range or original x range
    if(plot_common_x == TRUE) {
      modlist_plot <- modlist_valid_interp
    } else {
      modlist_plot <- modlist_valid
    }
    
    # Make Amoroso titles
    if(!is.null(amo_hell_cdf)) {
      amo_hell_cdf_title <- paste0("Amoroso"," (",amo_hell_cdf$method_short, ")")
    } else {
      amo_hell_cdf_title <- ""
    }
    if(!is.null(amo_hell_pdf)) {
      amo_hell_pdf_title <- paste0("Amoroso"," (",amo_hell_pdf$method_short, ")")  
    } else {
      amo_hell_pdf_title <- ""
    }
    
    # Make vector with all titles
    all_titles <- list(
      rdens = "R density()",
      scKDE_2infplus = "Adj. KDE ('twoInflections+')",
      mnorm = "Mixed Normal",
      amo_hell_cdf = amo_hell_cdf_title,
      amo_hell_pdf = amo_hell_pdf_title)
    
    # Extract only titles of valid models
    valid_titles <- all_titles[intersect(names(all_titles), names(modlist_valid))]

    # Make 3 of the titles more conscise
    if ("scKDE_2infplus" %in% names(valid_titles)) {
      valid_titles$scKDE_2infplus <- "Adj. KDE (2Inf+)"
    }
    
    for (key in c("amo_hell_cdf", "amo_hell_pdf")) {
      if (key %in% names(valid_titles)) {
        sign <- substr(valid_titles[[key]], nchar(valid_titles[[key]]) - 2, nchar(valid_titles[[key]]) - 2)
        valid_titles[[key]] <- sub("HELL-(CDF|PDF) \\(.*\\)", paste0("Hell-\\1", sign, ")"), valid_titles[[key]])
      }
    }
    
    # Tranform to vector
    titlevec <- valid_titles %>% unlist()
    
    # Create vector of colours
    colors <- c("royalblue", "springgreen4", "orange3","hotpink3","purple3")
    
    
    #---------------------------------------------------------------------------
    # Non-minimal style (i.e., not for proposal)
    #---------------------------------------------------------------------------
    
    if (minimal == FALSE) {
      
      # Initialize 1x5 plotting grid
      par(mfrow=c(1,5), oma = c(1, 1, 5, 1), cex.axis = 0.9, font.lab = 2,
          font.axis = 1, family = "Times New Roman")
      
      # Plot density estimates
      for (i in 1:length(modlist_plot)) {
        
        # Make empty plot
        plot(NA, xlim = c(xmin_plot, xmax_plot), ylim = c(0.0, ymax),
             xlab = "", ylab = "Density", main = titlevec[i], axes = FALSE)
        ifelse(is.null(xticks),
               axis(1),
               axis(1, at = xticks, labels = xticks))
        axis(2, las = 2)
        if(rug) {rug(dat, col = "blue", lwd = 1)}
        
        # Optional: add histogram
        if (hist == TRUE) {
          hist(dat, prob = T, breaks = breaks, col = "grey95",
               border = "grey85", axes = FALSE, add = TRUE)
        }
        
        # If model is valid, add density estimate
        if (length(modlist_plot[[i]]) > 1) {
          lines(modlist_plot[[i]]$x, modlist_plot[[i]]$y,
                col = colors[i], lwd = 1.7)
        }
        
        # Create empty label for data-generating distribution
        truedistlabel = NULL
        
        # Optional: add data-generating normal distribution
        if (length(generatingnormal==2) && is.numeric(generatingnormal)) {
          truedistlabel <- paste0("Normal(",
                                  "mean=", generatingnormal[1],", ",
                                  "sd=", generatingnormal[2],")")
          lines(xvals, dnorm(xvals,
                             mean = generatingnormal[1],
                             sd = generatingnormal[2]), 
                type = "l", lwd = 1, lty = 2, col = "grey30")
        }
        
        # Optional: add data-generating Amoroso distribution
        if (length(generatingamoroso==4) && is.numeric(generatingamoroso)) {
          truedistlabel <- paste0("Amoroso(",
                                  "a=", generatingamoroso[1],", ",
                                  "l=", generatingamoroso[2],", ",
                                  "c=", generatingamoroso[3],", ",
                                  "mu=", generatingamoroso[4],")")
          lines(xvals, dgg4(xvals,
                             a = generatingamoroso[1],
                             l = generatingamoroso[2],
                             c = generatingamoroso[3],
                             mu = generatingamoroso[4]),
                type = "l", lwd = 1, lty = 2, col = "grey30")
        }
      }
      
      # Add big title
      if (is.null(main)) {
        big_title <- "Nonparametric, Mixed Normal and Amoroso Fits"
        #amo_hell_pdf_title
      } else {
        big_title <- main
      }
      mtext(big_title, outer = TRUE, cex = 1.5, line = 2, font = 2)
      
      # Make function to add legend outside of the five plots
      add_legend <- function(...) {
        opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 1.2, 2), 
                    mar=c(0, 0, 0, 0), new=TRUE)
        on.exit(par(opar))
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend(...)
      }
      
      # Add true distribution legend
      if (!is.null(truedistlabel)) {
        add_legend("topright", legend=truedistlabel, cex=1.2,bty="n")
      }
      
      #---------------------------------------------------------------------------
      # Minimal style (for proposal)
      #---------------------------------------------------------------------------
      
    } else {
      
      # Initialize 1x5 plotting grid
      par(mfrow=c(1,5), mar = c(1,4,2,1), oma = c(2,0,2,0), cex.axis = 0.9, font.lab = 2,
          font.axis = 1, family = "Times New Roman", cex.main = 1.6)
      
      # Plot density estimates
      for (i in 1:length(modlist_plot)) {
        
        # Make empty plot
        plot(NA, xlim = c(xmin_plot, xmax_plot), ylim = c(0.0, ymax),
             xlab = "", ylab = "Density", main = titlevec[i], axes = FALSE)
        ifelse(is.null(xticks),
               axis(1),
               axis(1, at = xticks, labels = xticks))
        ifelse(is.null(yticks),
               axis(2, las = 2),
               axis(2, las = 2, at = yticks, labels = yticks))
            
        if(rug) {rug(dat, col = "cornflowerblue", lwd = 1)}
        
        # Optional: add histogram
        if (hist == TRUE) {
          hist(dat, prob = T, breaks = breaks, col = "grey95",
               border = "grey85", axes = FALSE, add = TRUE)
        }
        
        # If model is valid, add density estimate
        if (length(modlist_plot[[i]]) > 1) {
          lines(modlist_plot[[i]]$x, modlist_plot[[i]]$y,
                col = colors[i], lwd = 2)
        }
        
        # Create empty label for data-generating distribution
        truedistlabel = NULL
        
        # Optional: add data-generating normal distribution
        if (length(generatingnormal==2) && is.numeric(generatingnormal)) {
          truedistlabel <- paste0("Normal(",
                                  "mean=", generatingnormal[1],", ",
                                  "sd=", generatingnormal[2],")")
          lines(xvals, dnorm(xvals,
                             mean = generatingnormal[1],
                             sd = generatingnormal[2]), 
                type = "l", lwd = 1, lty = 2, col = "grey30")
        }
        
        # Optional: add data-generating Amoroso distribution
        if (length(generatingamoroso==4) && is.numeric(generatingamoroso)) {
          truedistlabel <- paste0("Amoroso(",
                                  "a=", generatingamoroso[1],", ",
                                  "l=", generatingamoroso[2],", ",
                                  "c=", generatingamoroso[3],", ",
                                  "mu=", generatingamoroso[4],")")
          lines(xvals, dgg4(xvals,
                            a = generatingamoroso[1],
                            l = generatingamoroso[2],
                            c = generatingamoroso[3],
                            mu = generatingamoroso[4]),
                type = "l", lwd = 1, lty = 2, col = "grey30")
        }
        
        # # Add big title at left side
        # if (is.null(main)) {
        #   big_title <- paste0("rnorm(", as.character(n), ")")
        # } else {
        #   big_title <- main
        # }
        # mtext(big_title, outer = TRUE, cex = 1.5, line = 2, font = 2)
      }
    }
    
  } else {
    
    print("no plot (since plot = FALSE")
    
  }
  
  #############################
  ### 4. RETURN MODEL LISTS ###
  #############################
  
  return_list <- list(modlist, modlist_valid, modlist_valid_interp)
  names(return_list) <- c("modlist", "modlist_valid", "modlist_valid_interp")
  invisible(return_list)
  
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### Test the function ###

#data <- palmerpenguins::penguins$bill_depth_mm
#dat <- palmerpenguins::penguins$bill_length_mm
#dat <- palmerpenguins::penguins$flipper_length_mm
#res <- estimate_amoroso_np(dat, hist = TRUE, minimal = FALSE)
#res$modlist_valid

#set.seed(125)
#data <- rgg4(30, a=4,l=1,c=7,mu=0)
#data <- rgg4(100, a=4,l=1,c=6,mu=0)
#data <- rnorm(70, mean = 4, sd = 0.7)
#data <- rgg4(20, a=4, l=1, c=7, mu=0)

#res1 <- estimate_methods(dat = data, plot_common_x = TRUE, generatingnormal = c(4,0.7))
#res2 <- estimate_methods(dat = data, plot_common_x = TRUE, generatingamoroso = c(4,1,7,0))
#res3 <- estimate_methods(dat = data, plot_common_x = TRUE)

#names(res1$modlist_valid_interp)
