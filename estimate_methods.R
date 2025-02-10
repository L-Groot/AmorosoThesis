#-------------------------------------------------------------------------------
# Load functions and packages
#-------------------------------------------------------------------------------

# -> for estimating Amoroso
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso.R"))
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
require(ggplot2)
require(tidyverse)
require(tidyr)
require(gridExtra)
require(gamlss.dist)


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
fit_amoroso <- function(dat, amoinaplus = FALSE) {
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

estimate_methods <- function(dat, amoinaplus = FALSE) {
  dat <- clean_data(dat)  # Remove NAs
  n <- length(dat)  # Get sample size
  
  # Estimate Amoroso
  amo <- fit_amoroso(dat, amoinaplus)
  
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
# plot_methods function
#-------------------------------------------------------------------------------

plot_methods <- function(dat, res,
                         plot_common_x = TRUE,
                         generatingnormal = NULL, #supply (mean,sd)
                         generatingamoroso = NULL, #supply (a,l,c,mu)
                         generatingexgauss = NULL, #supply (mu,sigma,tau)
                         xticks = NULL,
                         yticks = NULL,
                         ymax = NULL,
                         xmin = NULL,
                         xmax = NULL,
                         bins = 30) {
  
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
  
  cat("ymax: ", ymax, "\n")
  cat("xmin_plot: ", xmin_plot, "\n")
  cat("xmax_plot: ", xmax_plot, "\n")
  
  # Get titles and colors (from your function)
  all_titles <- list(
    rdens = "R density()",
    scKDE_2infplus = "Adj. KDE (2Inf+)",
    mnorm = "Mixed Normal",
    amo_hell_cdf = if (length(modlist_all$amo_hell_cdf) > 1) {
      sign <- ifelse(modlist_all$amo_hell_cdf$method_short == "HELL-CDF (+)", "+", "-")
      paste0("Amoroso (Hell-CDF", sign, ")")
    } else "",
    amo_hell_pdf = if (length(modlist_all$amo_hell_pdf) > 1) {
      sign <- ifelse(modlist_all$amo_hell_pdf$method_short == "HELL-CDF (+)", "+", "-")
      paste0("Amoroso (Hell-CDF", sign, ")")
    } else ""
  )
  
  # Extract valid model titles
  valid_titles <- all_titles[intersect(names(all_titles), names(modlist_valid))]
  
  # Convert to vector
  titlevec <- unlist(valid_titles)
  
  # Create vector of colours
  colors <- c("royalblue", "springgreen4", "orange3","hotpink3","purple3")
  
  # Convert modlist_plot into a tidy dataframe for ggplot
  df_list <- lapply(names(modlist_plot), function(name) {
    data.frame(
      x = modlist_plot[[name]]$x,
      y = modlist_plot[[name]]$y,
      method = name
    )
  })
  
  df_plot <- bind_rows(df_list)
  
  # Assign colors based on method names
  df_plot$color <- colors[match(df_plot$method, names(valid_titles))]
  
  # Create histogram data
  hist_data <- data.frame(x = dat)
  
  # Generate plots
  # Create a list to store the plots
  plot_list <- list()
  
  # Loop over each method in valid_titles
  for (meth in names(valid_titles)) {
    
    # Filter the data for the current method
    method_data <- df_plot %>% filter(method == meth)
    
    # Get color for current method
    color <- unique(method_data$color)
    
    # Create the plot for the current method
    p <- ggplot() +
      # Histogram
      geom_histogram(data = hist_data, aes(x = x, y = ..density..), 
                     bins = bins, fill = "grey95", color = "grey85", alpha = 0.5)
    
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
      labs(title = valid_titles[[meth]], x = NULL, y = "Density") +
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
  
  # Arrange in a 1-row, 5-column grid
  grid.arrange(grobs = plot_list, ncol = 5)
}





#-------------------------------------------------------------------------------
# Examples
#-------------------------------------------------------------------------------

# Define data
#dat <- palmerpenguins::penguins$flipper_length_mm
#dat <- palmerpenguins::penguins$bill_length_mm
#dat <- palmerpenguins::penguins$bill_depth_mm

             