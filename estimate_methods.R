#-------------------------------------------------------------------------------
# Load functions and packages
#-------------------------------------------------------------------------------

# -> for estimating Amoroso
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso_hell_aplus.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/mnorm_functions.R"))

# -> for Amoroso density function
require(AmoRosoDistrib)
# -> for adjusted KDE
require(scdensity) 
# -> for mixed normal estimation
require(mclust)
# -> for mixed normal density function
require(LaplacesDemon)

# Other packages
require(caret)
require(ggplot2)
require(tidyverse)
require(tidyr)
require(gridExtra)
require(gamlss.dist)
require(tidymodels)
require(spatstat)


#-------------------------------------------------------------------------------
# Helper functions
#-------------------------------------------------------------------------------

# Function to safely execute expressions and handle errors
safe_execute <- function(expr, object_name, data_vector) {
  tryCatch(
    eval(bquote(.(expr)), envir = list(dat = data_vector)),
    error = function(e) {
      cat(sprintf("Error with fitting %s: %s; Other methods were still fit.\n",
                  object_name, e$message))
      return(NA)  # Assigns NA if method fails
    }
  )
}


# Function to remove NAs and print a warning if necessary
clean_data <- function(dat) {
  num_nas_to_remove <- sum(is.na(dat))
  dat <- as.vector(na.omit(dat))
  
  if (num_nas_to_remove > 0) {
    message(sprintf("WARNING: %d NAs were removed from the data.\n", num_nas_to_remove))
    cat("--------------------------------------------------------------------\n")
  }
  return(dat)
}

# Function to estimate Amoroso distribution
fit_amoroso_all <- function(dat, amoinaplus = FALSE) {
  amo <- safe_execute(quote(estimate_amoroso(dat, plot = 0, criterion = "ML")), "amo", dat)
  
  if (length(amo)==1) return(NULL)  # Return NULL if estimation failed
  
  amo_x <- amo$x
  amo <- if (amoinaplus) {
    amo$all_models %>% filter(space == "+")  # Only a > 0 parameter space
  } else {
    amo$max_L_models  # Best likelihood model (a+ or a-)
  }
  
  return(list(amo = amo, amo_x = amo_x))
}

# Function to only estimate Hellinger Amorosos in a>0
fit_amoroso_hell_aplus <- function(dat) {
  amo <- safe_execute(quote(estimate_amoroso_hell_aplus(dat, plot = 0, criterion = "ML")), "amo", dat)
  
  if (length(amo)==1) return(NULL)  # Return NULL if estimation failed
  
  amo_x <- amo$x
  amo <- amo$max_L_models
  
  return(list(amo = amo, amo_x = amo_x))
}

# Function to estimate different density models
fit_others <- function(dat) {
  list(
    rdens = safe_execute(quote(density(dat)), "rdens", dat),
    scKDE_2infplus = safe_execute(quote(scdensity(dat, constraint = "twoInflections+")), "scKDE_2infplus", dat),
    mnorm = {
      mnorm_fit <- safe_execute(quote(densityMclust(dat, plot = FALSE)), "mnorm", dat)
      if (length(mnorm_fit)>1) {
        mnorm_fit$x <- density(dat)$x  # Use same x as base R density
        mnorm_fit$y <- predict_mnorm(mnorm_fit$x, mnorm_fit, plot = FALSE)
      }
      mnorm_fit
    }
  )
}

# Function to extract specific Amoroso fits
extract_amoroso <- function(amo, method_id) {
  amo_df <- amo$amo
  amo_x <- amo$amo_x
  if (!(method_id %in% amo_df$method_ID)) {
    warning(paste(method_id, "Amoroso absent"))
    return(NULL)
  }
  
  model <- amo_df %>% filter(method_ID == method_id)
  pars <- model %>% slice(1) %>% select(a, l, c, mu) %>% unlist(use.names = FALSE)
  density_values <- dgg4(amo_x, pars[1], pars[2], pars[3], pars[4])
  density_values[is.na(density_values)] <- 0
  
  return(list(x = amo_x, y = density_values, pars = pars,
              method = paste0(model$method, " (", model$space, ")"),
              method_short = paste0(model$method_ID, " (", model$space, ")")))
}

# Function to define plot limits based on density models
define_plot_limits <- function(dat, modlist_valid, breaks = 20) {
  xmin <- min(sapply(modlist_valid, function(mod) min(mod$x)))
  xmax <- max(sapply(modlist_valid, function(mod) max(mod$x)))
  xvals <- seq(xmin, xmax, length.out = 512)
  
  ymax <- max(sapply(modlist_valid, function(mod) max(mod$y)), na.rm = TRUE)
  
  # Make ymax (either prespecified or dynamically determined)
  hist_list <- hist(dat, breaks=breaks, plot=F)
  if(is.null(ymax)) {
    ymax <- max(hist_list$density)
  } else {
    ymax <- ymax
  }
  
  buffer <- 0.15*ymax
  ymax <- ymax + buffer
  
  # Define expanded x-range for plots
  xrange <- diff(range(dat))
  xmin_plot <- min(dat) - 0.2 * xrange
  xmax_plot <- max(dat) + 0.2 * xrange
  
  return(list(xvals = xvals, ymax = ymax, xmin_plot = xmin_plot, xmax_plot = xmax_plot))
}




#-------------------------------------------------------------------------------
# estimate_methods() function
#-------------------------------------------------------------------------------

estimate_methods <- function(dat) {
  dat <- clean_data(dat)  # Remove NAs
  n <- length(dat)  # Get sample size
  
  # Estimate Amoroso
  amo <- fit_amoroso_hell_aplus(dat)
  
  # Extract Hellinger models if Amoroso fit succeeded
  amo_hell_cdf <- extract_amoroso(amo, "HELL-CDF")
  amo_hell_pdf <- extract_amoroso(amo, "HELL-PDF")
  
  # Estimate density models
  other_three_models <- fit_others(dat)
  
  # Compile model lists
  modlist <- c(other_three_models, list(amo_hell_cdf = amo_hell_cdf, amo_hell_pdf = amo_hell_pdf))
  modlist_valid <- modlist[sapply(modlist, function(mod) length(mod) > 1)]
  
  # Define plot limits
  plot_limits <- define_plot_limits(dat, modlist_valid)
  xvals <- plot_limits$xvals
  ymax <- plot_limits$ymax
  xmin_plot <- plot_limits$xmin_plot
  xmax_plot <- plot_limits$xmax_plot
  
  # Interpolate valid models
  modlist_valid_interp <- lapply(modlist_valid, function(mod) {
    if (!is.null(mod)) {
      y <- approx(mod$x, mod$y, xout = xvals, rule = 1)$y
      y[is.na(y)] <- 0
      list(x = xvals, y = y)
    } else {
      NULL
    }
  })
  
  return(invisible(list(
    modlist = modlist,
    modlist_valid_interp = modlist_valid_interp,
    ymax = ymax,
    xmin_plot = xmin_plot,
    xmax_plot = xmax_plot
  )))
}




#-------------------------------------------------------------------------------
# plot_methods() function
#-------------------------------------------------------------------------------
plot_methods <- function(dat, res,
                         plot_common_x = TRUE,
                         generatingnormal = NULL, #supply (mean,sd)
                         generatingamoroso = NULL, #supply (a,l,c,mu)
                         generatingexgauss = NULL, #supply (mu,sigma,nu)
                         xticks = NULL,
                         yticks = NULL,
                         ymax = NULL,
                         xmin = NULL,
                         xmax = NULL,
                         bins = 30,
                         method_to_plot = NULL,
                         alpha = 1,
                         main = NULL) {  # New argument 'main'
  
  modlist_all <- res$modlist
  modlist_valid <- res$modlist_valid
  modlist_valid_interp <- res$modlist_valid_interp
  
  modlist_plot <- if (plot_common_x) modlist_valid_interp else modlist_valid
  
  ifelse(is.null(ymax), ymax <- res$ymax, ymax <- ymax)
  ifelse(is.null(xmin), xmin_plot <- res$xmin_plot, xmin_plot <- xmin)
  ifelse(is.null(xmax), xmax_plot <- res$xmax_plot, xmax_plot <- xmax)
  
  cat("ymax: ", ymax, "\n")
  cat("xmin_plot: ", xmin_plot, "\n")
  cat("xmax_plot: ", xmax_plot, "\n")
  
  all_titles <- list(
    rdens = "R density()",
    scKDE_2infplus = "Adj. KDE (2Inf+)",
    mnorm = "Mixed Normal",
    amo_hell_cdf = "Amoroso (Hell-CDF)",
    amo_hell_pdf = "Amoroso (Hell-PDF)"
  )
  
  valid_titles <- all_titles[intersect(names(all_titles), names(modlist_valid))]
  titlevec <- unlist(valid_titles)
  colors <- c("royalblue", "springgreen4", "orange3","hotpink3","purple3")
  
  df_list <- lapply(names(modlist_plot), function(name) {
    data.frame(
      x = modlist_plot[[name]]$x,
      y = modlist_plot[[name]]$y,
      method = name
    )
  })
  
  df_plot <- bind_rows(df_list)
  df_plot$color <- colors[match(df_plot$method, names(valid_titles))]
  
  hist_data <- data.frame(x = dat)
  
  if (!is.null(method_to_plot) && method_to_plot %in% names(valid_titles)) {
    method_data <- df_plot %>% filter(method == method_to_plot)
    color <- unique(method_data$color)
    
    p <- ggplot() +
      geom_histogram(data = hist_data, aes(x = x, y = ..density..),
                     bins = bins, fill = "grey90", color = "grey80", alpha = 0.5)
    
    if (!is.null(generatingnormal)) {
      p <- p + geom_line(aes(x = seq(xmin_plot, xmax_plot, length.out = 100),
                             y = dnorm(seq(xmin_plot, xmax_plot, length.out = 100),
                                       mean = generatingnormal[1], 
                                       sd = generatingnormal[2])), 
                         color = "grey60", linetype = "dashed", size = 1)
    }
    
    if (!is.null(generatingamoroso)) {
      p <- p + geom_line(aes(x = seq(xmin_plot, xmax_plot, length.out = 100),
                             y = dgg4(seq(xmin_plot, xmax_plot, length.out = 100),
                                      generatingamoroso[1], generatingamoroso[2], 
                                      generatingamoroso[3], generatingamoroso[4])), 
                         color = "grey60", linetype = "dashed", size = 1)
    }
    
    if (!is.null(generatingexgauss)) {
      p <- p + geom_line(aes(x = seq(xmin_plot, xmax_plot, length.out = 100),
                             y = dexGAUS(seq(xmin_plot, xmax_plot, length.out = 100),
                                         generatingexgauss[1], generatingexgauss[2], 
                                         generatingexgauss[3])), 
                         color = "grey60", linetype = "dashed", size = 0.7)
    }
    
    p <- p +
      geom_line(data = method_data, aes(x = x, y = y), color = color, size = 0.8, alpha = alpha) +
      labs(title = ifelse(is.null(main), valid_titles[[method_to_plot]], main),  # Use 'main' if provided
           x = NULL, y = "Density") +
      theme_minimal() +
      theme(
        text = element_text(family = "Times"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 12)),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")
      ) +
      coord_cartesian(xlim = c(xmin_plot, xmax_plot), ylim = c(0, ymax))
    
    if (!is.null(xticks)) {
      p <- p + scale_x_continuous(breaks = xticks)
    }
    
    if (!is.null(yticks)) {
      p <- p + scale_y_continuous(breaks = yticks)
    }
    
    print(p)
  } else {
    plot_list <- list()
    for (meth in names(valid_titles)) {
      method_data <- df_plot %>% filter(method == meth)
      color <- unique(method_data$color)
      p <- ggplot() +
        geom_histogram(data = hist_data, aes(x = x, y = ..density..),
                       bins = bins, fill = "grey90", color = "grey80", alpha = 0.5) +
        geom_line(data = method_data, aes(x = x, y = y), color = color, size = 0.8, alpha = alpha) +
        labs(title = ifelse(is.null(main), valid_titles[[meth]], main),  # Use 'main' if provided
             x = NULL, y = "Density") +
        theme_minimal() +
        theme(
          text = element_text(family = "Times"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 12)),
          axis.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 10),
          plot.margin = margin(1, 1, 1, 1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black")
        ) +
        coord_cartesian(xlim = c(xmin_plot, xmax_plot), ylim = c(0, ymax))
      plot_list[[meth]] <- p
    }
    grid.arrange(grobs = plot_list, ncol = 5)
  }
}



#-------------------------------------------------------------------------------
# Hell-CDF vs Hell-PDF function
#-------------------------------------------------------------------------------

hellcdf_vs_hellpdf <- function(dat) {
  
  # Get maxL Amorosos (for each method either +ve or_ve parameter space)
  amo_maxL_df <- estimate_amoroso_hell_aplus(dat)$max_L_models
  
  # Identify whether Amoroso Hell-CDF or Hell-PDF has higher likelihood
  hell_pdf_row <- which(amo_maxL_df$method_ID == "HELL-PDF")
  hell_cdf_row <- which(amo_maxL_df$method_ID == "HELL-CDF")
  
  # Compare their positions
  if (hell_pdf_row < hell_cdf_row) {
    # If Hell-PDF Amo has higher likelihood, choose that one
    maxL_amo <- "amo_hell_pdf"
    sign <- amo_maxL_df$space[hell_pdf_row]
    amo_title <- paste0("Amoroso (Hell-PDF", sign, ")")
    amo_color <- "hotpink3"
    print("HELL-PDF has higher maxL than HELL-CDF.")
    return(list(maxL_amo = maxL_amo, sign = sign, amo_title = amo_title,
                amo_color = amo_color))
    
  } else if (hell_pdf_row > hell_cdf_row) {
    # If Hell-CDF Amo has higher likelihood, choose that one
    maxL_amo <- "amo_hell_cdf"
    sign <- amo_maxL_df$space[hell_cdf_row]
    amo_title <- paste0("Amoroso (Hell-CDF", sign, ")")
    amo_color <- "purple3"
    print("HELL-CDF appears has higher maxL than HELL-PDF.")
    return(list(maxL_amo = maxL_amo, sign = sign, amo_title = amo_title,
                amo_color = amo_color))
    
  } else {
    # Else throw error cuz something must have gone wrong:(
    stop("Error when comparing hell-cdf and hell-pdf row positions")
  }
  
}




#-------------------------------------------------------------------------------
# plot_rdens_amo() function
#-------------------------------------------------------------------------------

plot_rdens_amo <- function(dat, res,
                         plot_common_x = TRUE,
                         generatingnormal = NULL, #supply (mean,sd)
                         generatingamoroso = NULL, #supply (a,l,c,mu)
                         generatingexgauss = NULL, #supply (mu,sigma,nu)
                         xticks = NULL,
                         yticks = NULL,
                         ymax = NULL,
                         xmin = NULL,
                         xmax = NULL,
                         rug = T,
                         bins = 30) {
  
  # dat <- exGauss_simdat
  # xmin <- 130
  # xmax <- 530
  # ymax <- 0.01
  
  # Get valid models
  modlist_all <- res$modlist
  modlist_valid <- res$modlist_valid
  modlist_valid_interp <- res$modlist_valid_interp
  
  # Choose whether to interpolate x-range across models
  modlist_plot <- if (plot_common_x) modlist_valid_interp else modlist_valid
  
  # Determine plot bounds
  ifelse(is.null(ymax),
         ymax <- res$ymax,
         ymax <- ymax)
  ifelse(is.null(xmin),
         xmin_plot <- res$xmin_plot,
         xmin_plot <- xmin)
  ifelse(is.null(xmax),
         xmax_plot <- res$xmax_plot,
         xmax_plot <- xmax)
  
  # Identify whether Hell-CDF or Hell-PDF has higher likelihood for data
  maxL_amo <- hellcdf_vs_hellpdf(dat)
  
  # Make vector of titles
  plot_titles <- c("R density()", maxL_amo$amo_title)
  
  # Extract only R density and maxL Hellinger Amoroso
  modlist_plot <- res$modlist_valid_interp[c("rdens", maxL_amo$maxL_amo)]
  
  # Convert modlist to df for ggplot
  df_list <- lapply(names(modlist_plot), function(name) {
    data.frame(
      x = modlist_plot[[name]]$x,
      y = modlist_plot[[name]]$y,
      method = name
    )
  })
  
  df_plot <- bind_rows(df_list)
  
  # Specify two colours
  two_colors <- c("royalblue", maxL_amo$amo_color)
  
  # Assign colors based on method names
  df_plot$color <- two_colors[match(df_plot$method, names(modlist_plot))]
  
  # Make the two fit plots

  # Create histogram data
  hist_data <- data.frame(x = dat)
  
  # Create a list to store the plots
  plot_list <- list()
  
  
  # Loop over each method in valid_titles
  for (i in 1:2) {
    
    # Filter the data for the current method
    meth <- names(modlist_plot)[i]
    method_data <- df_plot %>% filter(method == meth)
    
    # Get color for current method
    color <- unique(method_data$color)
    
    # Create the plot for the current method
    p <- ggplot() +
      # Histogram
      geom_histogram(data = hist_data, aes(x = x, y = ..density..), 
                     bins = bins, fill = "grey95", color = "grey85", alpha = 0.5)
    # Optional: Add rug marks to show individual data points
    if (rug) {
      p <- p + geom_rug(data = hist_data, aes(x = x), color = "black", alpha = 0.6)
    }
    
    # Add generatingnormal line if provided
    if (!is.null(generatingnormal)) {
      p <- p + geom_line(aes(x = seq(xmin_plot, xmax_plot, length.out = 100), 
                             y = dnorm(seq(xmin_plot, xmax_plot, length.out = 100), 
                                       mean = generatingnormal[1], 
                                       sd = generatingnormal[2])), 
                         color = "grey", linetype = "dashed", size = 1)
    }
    
    # Add generatingamoroso line if provided
    if (!is.null(generatingamoroso)) {
      # Assuming you have the amoroso distribution function, let's call it `dgg4()`
      p <- p + geom_line(aes(x = seq(xmin_plot, xmax_plot, length.out = 100), 
                             y = dgg4(seq(xmin_plot, xmax_plot, length.out = 100), 
                                      generatingamoroso[1], generatingamoroso[2], 
                                      generatingamoroso[3], generatingamoroso[4])), 
                         color = "grey", linetype = "dashed", size = 1)
    }
    
    # Add generatingexgauss line if provided
    if (!is.null(generatingexgauss)) {
      # Assuming you have the amoroso distribution function, let's call it `dgg4()`
      p <- p + geom_line(aes(x = seq(xmin_plot, xmax_plot, length.out = 100), 
                             y = dexGAUS(seq(xmin_plot, xmax_plot, length.out = 100), 
                                         generatingexgauss[1], generatingexgauss[2], 
                                         generatingexgauss[3])), 
                         color = "grey70", linetype = "dashed", size = 0.7)
    }
    
    # Add density estimate
    p <- p + 
      geom_line(data = method_data, 
                aes(x = x, y = y), color = color, size = 0.8) +
      
      # Add other stuff
      labs(title = plot_titles[[i]], x = NULL, y = "Density") +
      theme_minimal() +
      theme(
        text = element_text(family = "Times"),  # Use Times New Roman
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 12)),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")
      ) +
      coord_cartesian(xlim = c(xmin_plot, xmax_plot), ylim = c(0, ymax))
    
    # Apply xticks if not NULL
    if (!is.null(xticks)) {
      p <- p + scale_x_continuous(breaks = xticks)
    }
    
    # Apply yticks if not NULL
    if (!is.null(yticks)) {
      p <- p + scale_y_continuous(breaks = yticks)
    }
    
    # Add the plot to the list
    plot_list[[meth]] <- p
    
  }
  
  # Arrange in a 1-row, 2-column grid
  grid.arrange(grobs = plot_list, ncol = 2)
  
  # Return the plot list
  return(plot_list)
  
}





#-------------------------------------------------------------------------------
# plot_qq() function
#-------------------------------------------------------------------------------

theoretical_qq <- function(dat, res, method_id = "rdens",
                    generatingnormal = NULL,
                    generatingamoroso = NULL,
                    generatingexgauss = NULL,
                    lower_tail = TRUE,
                    par_on_axislab = FALSE,
                    rev = TRUE) {
  
  # dat = exGauss_simdat
  # res = estimate_methods(dat)
  # method_id = "amo_hell_pdf"
  # generatingnormal = NULL
  # generatingamoroso = NULL
  # generatingexgauss = c(mu,sigma,nu)
  # lower_tail = TRUE
  # par_on_axislab = FALSE
  # rev = F
  
  # Check if exactly one of the generating arguments is not NULL
  non_null_args <- sum(!sapply(list(generatingnormal, generatingamoroso,
                                    generatingexgauss), is.null))
  
  if (non_null_args != 1) {
    stop(paste0("Exactly one of generatingnormal, generatingamoroso, or",
                "generatingexgauss must be provided, but not more than one."))
  }
  
  #Extract target fit
  modlist <- res$modlist
  fit <- modlist[[method_id]]
  names(fit)
  
  # Get fit density values
  df <- data.frame(x = fit$y)
  
  # Get fit quantiles
  if(method_id == "rdens") {
    fit_quantiles <- quantile(density(dat), probs = seq(0.01,0.99,length.out=50))
    ylab <- "Fit quantiles: R density()"
  } else if (startsWith(method_id, "amo")) {
    pars <- fit$pars
    fit_quantiles <- qgg4(seq(0.01,0.99,length.out=50), pars[1],pars[2],pars[3],pars[4],
                          lower.tail = lower_tail)
    if(par_on_axislab) {
      ylab <- paste0("Fit quantiles: Amoroso(", round(pars[1],1),", ",
                     round(pars[2],1),", ", round(pars[3],1),", ",
                     round(pars[4],1),")")
    } else {
      ylab <- "Fit quantiles: Amoroso"
    }
  }
  if(rev) fit_quantiles <- rev(fit_quantiles)
  
  # Get theoretical quantiles (from data-generating distribution)
  if (!is.null(generatingnormal)) {
    pars <- generatingnormal
    true_quantiles <- qnorm(seq(0.01,0.99,length.out=50), pars[1], pars[2])
    if (par_on_axislab) {
      xlab <- paste0("True quantiles: Normal(", round(pars[1],1),", ",
                     round(pars[2],1), ")")
    } else {
      xlab <- "True quantiles: Normal"
    }
    
  } else if (!is.null(generatingamoroso)) {
    pars <- generatingamoroso
    true_quantiles <- qgg4(seq(0.01,0.99,length.out=50), pars[1], pars[2], pars[3], pars[4])
    if (par_on_axislab) {
      xlab <- paste0("True quantiles: Amoroso(", round(pars[1],1),", ",
                     round(pars[2],1),", ", round(pars[3],1),", ",
                     round(pars[4],1),")")
    } else {
      xlab <- "True quantiles: Amoroso"
    }
    
  } else if (!is.null(generatingexgauss)) {
    pars <- generatingexgauss
    true_quantiles <- qexGAUS(seq(0.01,0.99,length.out=50), pars[1], pars[2], pars[3])
    if (par_on_axislab) {
      xlab <- paste0("True quantiles: ExGauss(", round(pars[1],1),", ",
                     round(pars[2],1),", ", round(pars[3],1), ")")
    } else {
      xlab <- "True quantiles: ExGauss"
    }
    
  }
  
  # Function to replace Inf/-Inf with 0 and throw a warning
  replace_inf_with_warning <- function(vec, name) {
    if (any(!is.finite(vec))) {
      warning(paste("Warning:", name, "contains Inf or -Inf. Replacing with 0."))
      vec[!is.finite(vec)] <- 0
    }
    return(vec)
  }
  
  # Replace any Inf and -Inf in true_quantiles and fit_quantiles with 0
  (true_quantiles <- replace_inf_with_warning(true_quantiles, "true_quantiles"))
  (fit_quantiles <- replace_inf_with_warning(fit_quantiles, "fit_quantiles"))
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Fit_Quantiles = fit_quantiles,
    True_Quantiles = true_quantiles
  )
  
  # Make Q-Q plot
  p <- ggplot(plot_data, aes(x = True_Quantiles, y = Fit_Quantiles)) +
          geom_point(size = 0.5, color = "chartreuse4") +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkorange3") +
          labs(
            title = "Q-Q plot",
            x = xlab,
            y = ylab
          ) +
          theme(
            text = element_text(family = "Times"),  # Use Times New Roman
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 12)),
            axis.title = element_text(size = 10, face = "bold"),
            axis.text = element_text(size = 10),
            plot.margin = margin(1, 1, 1, 1, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black")
          ) +
          coord_obs_pred() # same scale for x and y with aspect ratio 1
  
  # Return plot
  return(p)
}

#-------------------------------------------------------------------------------
# plot_pp() function
#-------------------------------------------------------------------------------

theoretical_pp <- function(dat, res, method_id = "rdens",
                           generatingnormal = NULL,
                           generatingamoroso = NULL,
                           generatingexgauss = NULL,
                           lower_tail = TRUE,
                           par_on_axislab = FALSE) {
  
  # method_id = "amo_hell_pdf"
  # generatingnormal = NULL
  # generatingamoroso = NULL
  # generatingexgauss = c(mu,sigma,nu)
  # lower_tail = TRUE
  # par_on_axislab = FALSE
  
  # Check if exactly one of the generating arguments is not NULL
  non_null_args <- sum(!sapply(list(generatingnormal, generatingamoroso,
                                    generatingexgauss), is.null))
  
  if (non_null_args != 1) {
    stop(paste0("Exactly one of generatingnormal, generatingamoroso, or",
                "generatingexgauss must be provided, but not more than one."))
  }
  
  # Extract target fit
  modlist <- res$modlist
  fit <- modlist[[method_id]]
  
  # Get fit density values
  df <- data.frame(x = fit$y)
  
  # Define x range
  #x <- seq(min(exGauss_simdat), max(exGauss_simdat), length.out = 50)
  x_min <- quantile(exGauss_simdat, 0.10)
  x_max <- quantile(exGauss_simdat, 0.90)
  
  # Generate sequence between the 10th and 90th percentiles
  x <- seq(x_min, x_max, length.out = 50)
  
  ### FIT CDF probabilites ###
  
  if(method_id == "rdens") {
    print("Fit is R density()")
    # Approximate function using linear interpolation
    f <- approxfun(fit$x, fit$y, method = "linear", yleft = 0, yright = 0)
    # Get integral (CDF prob) at each x value
    fit_percentiles <- c()
    for (xval in x) {
      integral <- integrate(f, -Inf, xval)$value
      fit_percentiles <- append(fit_percentiles, integral)
    }
    # Replace zeroes in tail with 1
    # Find the position of the last nonzero element
    #last_nonzero <- max(which(fit_percentiles != 0))
    # Replace all trailing zeros with 1
    #fit_percentiles[(last_nonzero + 1):length(fit_percentiles)] <- 1
    ylab <- "Fit percentile: R density()"
    
  } else if (startsWith(method_id, "amo")) {
    print("Fit is Amoroso")
    # Extract Amoroso fit parameters
    pars <- fit$pars
    # Get integral (CDF prob) at each x value
    fit_percentiles <- AmoRosoDistrib::pgg4(x, pars[1], pars[2], pars[3], pars[4])
    if (par_on_axislab) {
      ylab <- paste0("Fit percentiles: Amoroso(", round(pars[1],1),", ",
                     round(pars[2],1),", ", round(pars[3],1),", ",
                     round(pars[4],1),")")
    } else {
      ylab <- "Fit percentiles: Amoroso"
    } 
  }
  
  ### TRUE CDF probabilites ###
  
  # True Normal
  if (!is.null(generatingnormal)) {
    print("True is Normal")
    pars <- c(generatingnormal[1], generatingnormal[2])
    true_percentiles <- pnorm(x, pars[1], pars[2], lower.tail = lower_tail)
    if (par_on_axislab) {
      xlab <- paste0("True percentiles: Normal(", round(pars[1],1),", ",
                     round(pars[2],1), ")")
    } else {
      xlab <- "True percentiles: Normal"
    }
    
  
  # True Amoroso
  } else if (!is.null(generatingamoroso)) {
    print("True is Amoroso")
    pars <- c(generatingamoroso[1], generatingamoroso[2], generatingamoroso[3],
              generatingamoroso[4])
    true_percentiles <- AmoRosoDistrib::pgg4(x, pars[1], pars[2], pars[3], pars[4],
                                          lower.tail = lower_tail)
    if (par_on_axislab) {
      xlab <- paste0("True percentiles: Amoroso(", round(pars[1],1),", ",
                     round(pars[2],1),", ", round(pars[3],1),", ",
                     round(pars[4],1),")")
    } else {
      xlab <- "True percentiles: Amoroso"
    }
    
    
  # True exGauss
  } else if (!is.null(generatingexgauss)) {
    print("True is exGauss")
    pars <- c(generatingexgauss[1], generatingexgauss[2], generatingexgauss[3])
    true_percentiles <- gamlss.dist::pexGAUS(x, pars[1], pars[2], pars[3],
                                             lower.tail = lower_tail)
    if (par_on_axislab) {
      xlab <- paste0("True percentiles: ExGauss(", round(pars[1],1),", ",
                     round(pars[2],1),", ", round(pars[3],1), ")")
    } else {
      xlab <- "True percentiles: ExGauss"
    }
  }
  
  ### Make P-P plot ###
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Fit_Percentiles = fit_percentiles,
    True_Percentiles = true_percentiles
  )

  # Plot
  p <- ggplot(plot_data, aes(x = True_Percentiles, y = Fit_Percentiles)) +
    geom_point(size = 0.5, color = "chartreuse4") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkorange3") +
    labs(
      title = "P-P plot",
      x = xlab,
      y = ylab
    ) +
    theme(
      text = element_text(family = "Times"),  # Use Times New Roman
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 12)),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      plot.margin = margin(1, 1, 1, 1, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    coord_obs_pred() # same scale for x and y with aspect ratio 1
  
  # Return plot
  return(p)

}

#-------------------------------------------------------------------------------
# Examples
#-------------------------------------------------------------------------------

# xmin <- 130
# xmax <- 530
# ymax <- 0.01

# Define data
#dat <- palmerpenguins::penguins$flipper_length_mm
#dat <- palmerpenguins::penguins$bill_length_mm
#dat <- palmerpenguins::penguins$bill_depth_mm

# dat <- rnorm(40)
# res <- estimate_methods(dat)
# #plot_methods(dat, res)
# #plot_rdens_amo(dat, res, generatingnormal = c(0,1))
# theoretical_qq(dat, res, method_id = "rdens", generatingnormal = c(0,1))
# theoretical_qq(dat, res, method_id = "amo_hell_pdf", generatingnormal = c(0,1), lower_tail = T)
# theoretical_qq(dat, res, method_id = "amo_hell_cdf", generatingnormal = c(0,1), lower_tail = TRUE)

