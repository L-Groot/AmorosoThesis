################################################################################
# SIMULATIONS THAT INVESTIGATE THE PREDICTIVE PERFORMANCE OF THE FIVE METHODS  #
################################################################################

### Source functions ###
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/get_pp.R"))

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Helper function that plots one density
plot_1_density <- function(xvec, yvec, xmin, xmax, ymax, main = "", line_color = "steelblue",
                           yticks = NULL, xticks = NULL, lwd = 2, cex.tick = 1, yaxis = TRUE) {
  
  # Prepare data for ggplot2
  density_data <- data.frame(
    x = xvec,
    y = yvec
  )
  
  if (is.null(ymax)) {
    ymax <- max(density_data$y) * 1.1
  }
  
  p <- ggplot(density_data, aes(x = x, y = y)) +
    geom_line(size = lwd, color = line_color) +  # Custom line color and thickness
    
    # Styling
    labs(title = main, x = NULL, y = "Density") +
    
    # Theme and other customizations
    theme_minimal() +
    theme(
      text = element_text(family = "Times"),  # Times New Roman font
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 12)),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10 * cex.tick),  # Adjust tick font size based on cex.tick
      plot.margin = margin(1, 1, 1, 1, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      
      # Remove the legend
      legend.position = "none"
    )
  
  # If yticks is provided, modify the y-axis ticks
  if (!is.null(yticks)) {
    breaks <- seq(0, ymax, length.out = length(yticks))
    p <- p + scale_y_continuous(limits = c(0, ymax),
                                breaks = breaks,
                                labels = yticks)
  }
  
  # If xticks is provided, modify the x-axis ticks
  if (!is.null(xticks)) {
    breaks <- seq(xmin, xmax, length.out = length(xticks))
    p <- p + scale_x_continuous(limits = c(xmin, xmax),
                                breaks = breaks,
                                labels = xticks)
  }
  
  # Optionally remove the y-axis
  if (!yaxis) {
    p <- p + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank())
  }
  
  # Show the plot
  print(p)
  
  # Return the plot
  return(p)
}


################################################################################
### Define and plot distributions
################################################################################
#--------------
### Normal ###
#--------------

# Define parameters for the normal distributions
par_norm_wide <- c(100, 15)
par_norm_medium <- c(100, 10)
par_norm_narrow <- c(100, 5)

# Generate x values
x_vals <- seq(50, 150, length.out = 1000)

# Calculate normal density values
dens_wide <- dnorm(x_vals, mean = par_norm_wide[1], sd = par_norm_wide[2])
dens_medium <- dnorm(x_vals, mean = par_norm_medium[1], sd = par_norm_medium[2])
dens_narrow <- dnorm(x_vals, mean = par_norm_narrow[1], sd = par_norm_narrow[2])

# Create main titles for plots
lab1 <- expression(bold(paste("Normal(",
                              bolditalic("μ"), " = ", "100", ", ",
                              bolditalic("sd"), " = ", "15", ")")))
lab2 <- expression(bold(paste("Normal(",
                              bolditalic("μ"), " = ", "100", ", ",
                              bolditalic("sd"), " = ", "10", ")")))
lab3 <- expression(bold(paste("Normal(",
                              bolditalic("μ"), " = ", "100", ", ",
                              bolditalic("sd"), " = ", "5", ")")))

# Prepare parameters for plots
xmin = 50
xmax = 150
ymax = 0.08
yticks = c(0,0.08)
cex.tick = 2.5

# Plot the three Normals in separate plots
plot_1_density(x_vals, dens_wide, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick, yticks = c(0,0.08),
               line_color = "steelblue3")
plot_1_density(x_vals, dens_medium, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick, yticks = c(0,0.08),
               line_color = "steelblue3")
plot_1_density(x_vals, dens_narrow, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick, yticks = c(0,0.08),
               line_color = "steelblue3")

#--------------
### Amoroso ###
#--------------

# Define parameters for three Amoroso distributions
par_amo1 <- c(a = 2, l = 0.3, c = 5, mu = 0)
par_amo2 <- c(a = 5, l = 0.3, c = 5, mu = 0)
par_amo3 <- c(a = 7, l = 1.3, c = 8, mu = 0)

# Generate x values for the plot
x_vals <- seq(0, 9, length.out = 1000)

# Calculate Amoroso density values
dens_amo1 <- dgg4(x_vals, par_amo1[1], par_amo1[2], par_amo1[3], par_amo1[4])
dens_amo2 <- dgg4(x_vals, par_amo2[1], par_amo2[2], par_amo2[3], par_amo2[4])
dens_amo3 <- dgg4(x_vals, par_amo3[1], par_amo3[2], par_amo3[3], par_amo3[4])

# Make main title for each plot
lab1 <- expression(bold(paste("Amoroso(",
                         bolditalic("α"), " = ", "2", ", ", 
                         bolditalic("l"), " = ", "0.3", ", ", 
                         bolditalic("c"), " = ", "5", ", ", 
                         bolditalic("μ"), " = ", "0",
                         ")")))

lab2 <- expression(bold(paste("Amoroso(",
                         bolditalic("α"), " = ", "5", ", ", 
                         bolditalic("l"), " = ", "0.3", ", ", 
                         bolditalic("c"), " = ", "5", ", ", 
                         bolditalic("μ"), " = ", "0"),
                         ")"))

lab3 <- expression(bold(paste("Amoroso(",
                         bolditalic("α"), " = ", "7", ", ", 
                         bolditalic("l"), " = ", "1.3", ", ", 
                         bolditalic("c"), " = ", "8", ", ", 
                         bolditalic("μ"), " = ", "0",
                         ")")))

# Prepare parameters for plots
xmin = 0
xmax = 9
ymax = 0.65
yticks = c(0, 0.6)
xticks = c(0,3,6,9)
cex.tick = 2.5

# Plot the three Amorosos in separate plots
plot_1_density(x_vals, dens_amo1, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "springgreen3")
plot_1_density(x_vals, dens_amo2, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "springgreen3")
plot_1_density(x_vals, dens_amo3, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "springgreen3")

#------------------
### Ex-Gaussian ###
#------------------

# Define parameters for three ex-Gauss distributions
par_exgauss1 <- c(1, 0.1, 2)
par_exgauss2 <- c(1, 0.5, 2)
par_exgauss3 <- c(1, 0.9, 3)

# Generate x values for the plot
x_vals <- seq(-1, 10, length.out = 1000)

# Generate ex-Gaussian density values
dens_exgauss1 <- dexGAUS(x_vals, par_exgauss1[1], par_exgauss1[2],
                         par_exgauss1[3])
dens_exgauss2 <- dexGAUS(x_vals, par_exgauss2[1], par_exgauss2[2],
                         par_exgauss2[3])
dens_exgauss3 <- dexGAUS(x_vals, par_exgauss3[1], par_exgauss3[2],
                         par_exgauss3[3])

# Make main titles for plots
lab1 <- expression(bold(paste("ex-Gaussian(",
                              bolditalic("μ"), " = ", "1", ", ", 
                              bolditalic("σ"), " = ", "0.1", ", ", 
                              bolditalic("τ"), " = ", "2", 
                              ")")))
lab2 <- expression(bold(paste("ex-Gaussian(",
                              bolditalic("μ"), " = ", "1", ", ", 
                              bolditalic("σ"), " = ", "0.5", ", ", 
                              bolditalic("τ"), " = ", "2", 
                              ")")))
lab3 <- expression(bold(paste("ex-Gaussian(",
                              bolditalic("μ"), " = ", "1", ", ", 
                              bolditalic("σ"), " = ", "0.9", ", ", 
                              bolditalic("τ"), " = ", "3", 
                              ")")))

# Prepare parameters for plots
xmin = 0
xmax = 10
ymax = 0.5
yticks = c(0, 0.5)
xticks = c(0, 2.5, 5, 7.5, 10)
cex.tick = 2.5

# Plot the three ex-Gaussians in separate plots
plot_1_density(x_vals, dens_exgauss1, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "violetred3")
plot_1_density(x_vals, dens_exgauss2, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "violetred3")
plot_1_density(x_vals, dens_exgauss3, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "violetred3")

#-------------------
### Mixed Normal ###
#-------------------

# Define parameters for three mixed Normal distributions
par_mnorm1 <- c(0.6, 3, 1.2, 7, 1)
par_mnorm2 <- c(0.5, 3.7, 0.7, 6, 0.8)
par_mnorm3 <- c(0.8, 3.5, 1.5, 7, 1)

# Generate x values for the plot
x_vals <- seq(-1, 10, length.out = 1000) #0.9

# Generate mixed normal density values
dens_mnorm1 <- dmixnorm(x_vals, p=par_mnorm1[1],
                        mu1=par_mnorm1[2], sd1=par_mnorm1[3],
                        mu2=par_mnorm1[4], sd2=par_mnorm1[5])
dens_mnorm2 <- dmixnorm(x_vals, p=par_mnorm2[1],
                        mu1=par_mnorm2[2], sd1=par_mnorm2[3],
                        mu2=par_mnorm2[4], sd2=par_mnorm2[5])
dens_mnorm3 <- dmixnorm(x_vals, p=par_mnorm3[1],
                        mu1=par_mnorm3[2], sd1=par_mnorm3[3],
                        mu2=par_mnorm3[4], sd2=par_mnorm3[5])

# Make main titles for plots
lab1 <- expression(bold(paste("ex-Gaussian(",
                              bolditalic("p") [bold("1")], " = ", "0.6", ", ",
                              bolditalic("μ") [bold("1")], " = ", "3", ", ", 
                              bolditalic("sd") [bold("1")], " = ", "1.2", ", ", 
                              bolditalic("μ") [bold("2")], " = ", "7", ", ",
                              bolditalic("sd") [bold("2")], " = ", "1", 
                              ")")))

lab2 <- expression(bold(paste("ex-Gaussian(",
                              bolditalic("p") [bold("1")], " = ", "0.5", ", ",
                              bolditalic("μ") [bold("1")], " = ", "3.7", ", ", 
                              bolditalic("sd") [bold("1")], " = ", "0.7", ", ", 
                              bolditalic("μ") [bold("2")], " = ", "6", ", ",
                              bolditalic("sd") [bold("2")], " = ", "0.8", 
                              ")")))

lab3 <- expression(bold(paste("ex-Gaussian(",
                              bolditalic("p") [bold("1")], " = ", "0.8", ", ",
                              bolditalic("μ") [bold("1")], " = ", "3.5", ", ", 
                              bolditalic("sd") [bold("1")], " = ", "1.5", ", ", 
                              bolditalic("μ") [bold("2")], " = ", "7", ", ",
                              bolditalic("sd") [bold("2")], " = ", "1", 
                              ")")))



# Prepare parameters for plots
xmin = -1
xmax = 10
ymax = 0.3
yticks = c(0, 0.3)
xticks = c(0,2.5,5,7.5,10)
cex.tick = 2.5

# Plot the three mixed Normals in separate plots
plot_1_density(x_vals, dens_mnorm1, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "darkorange2")
plot_1_density(x_vals, dens_mnorm2, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "darkorange2")
plot_1_density(x_vals, dens_mnorm3, xmin = xmin, xmax = xmax, ymax = ymax,
               yaxis = F, main = "", cex.tick = cex.tick,
               yticks = yticks, xticks = xticks,
               line_color = "darkorange2")





################################################################################
### Analyze results: Function to make result dataframes
################################################################################

table(res_exgauss$pars1$n200$win_df$max_logL_method)

res <- res_exgauss
pars <- "pars1"

# Function to compute proportions and create a data frame
compute_proportions <- function(res, pars) {
  sample_sizes <- c("n25", "n50", "n100", "n200")
  
  # Create a list to store proportion tables
  prop_list <- lapply(sample_sizes, function(n) {
    table(res[[pars]][[n]]$win_df$max_logL_method) / 100
  })
  
  # Fill in any methods with 0 entries
  prop_df <- combine_prop_list(prop_list)
  
  colnames(prop_df) <- sample_sizes
  
  return(prop_df)
}


# Function to standardize row names and convert to a data frame
combine_prop_list <- function(prop_list) {
  # Define the required row names (method names in the correct order)
  all_methods <- c("amo_hell_cdf", "amo_hell_pdf", "mnorm", "rdens", "scKDE_2infplus")
  
  # Iterate over each list element, ensuring all methods are present and ordered
  prop_list_filled <- lapply(prop_list, function(x) {
    x <- setNames(as.numeric(x), names(x))  # Ensure numeric values
    x[setdiff(all_methods, names(x))] <- 0  # Add missing methods with 0
    x[all_methods]  # Reorder to match the desired order
  })
  
  # Combine into a data frame
  prop_df <- as.data.frame(do.call(cbind, prop_list_filled))
  
  # Set row names
  rownames(prop_df) <- all_methods
  
  return(prop_df)
}



################################################################################
### Analyze results: Normal distribution
################################################################################

res_normal <- readRDS("res_normal_vipasha.rds")

(norm1_df <- compute_proportions(res_normal, "pars1"))
(norm2_df <- compute_proportions(res_normal, "pars2"))
(norm3_df <- compute_proportions(res_normal, "pars3"))


################################################################################
### Analyze results: ex-Gaussian distribution
################################################################################

res_exgauss <- readRDS("res_exgauss.rds")

(exgauss1_df <- compute_proportions(res_exgauss, "pars1"))
(exgauss2_df <- compute_proportions(res_exgauss, "pars2"))
(exgauss3_df <- compute_proportions(res_exgauss, "pars3"))


################################################################################
### Analyze results: Amoroso distribution
################################################################################

res_amo <- readRDS("res_amo.rds")

(amo1_df <- compute_proportions(res_amo, "pars1"))
(amo2_df <- compute_proportions(res_amo, "pars2"))
(amo3_df <- compute_proportions(res_amo, "pars3"))


################################################################################
### Analyze results: Mixed normal distribution
################################################################################

res_mnorm <- readRDS("res_mnorm.rds")

(mnorm1_df <- compute_proportions(res_mnorm, "pars1"))
(mnorm2_df <- compute_proportions(res_mnorm, "pars2"))
(mnorm3_df <- compute_proportions(res_mnorm, "pars3"))

  
