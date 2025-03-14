# Set working directory to script location (only works in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required libraries
library(magick)
library(dplyr)
library(ggplot2)
library(rsvg)

# Source external function
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

# # Define before
# newfilename <- "mygif"
# batch_size <- 10
# max_y <- 0.2
# svgwidth <- 17
# svgheight <- 3.4
# pngwidth <- 3092

# Make data
# Read and preprocess data
# s01gs <- read.table("noisedat/S01GS.DAT") %>%
#   rename(
#     key = V1,
#     rt = V2,
#     ignore = V3,
#     odd = V4,
#     rsi = V5,
#     stim = V6
#   ) %>%
#   filter(rt >= 500 & rt <= 2500)
# 
# dat<- s01gs$rt
# 
# # Shuffle data for random sampling
# set.seed(123)  # For reproducibility
# dat <- sample(dat)[1:30]


#-------------------------------------------------------------------------------
make_gif <- function(dat = NULL,
                     newfilename = "test",
                     batch_size = 10,
                     max_y = NULL,
                     svgwidth = 17,
                     svgheight = 5,
                     pngwidth = 3092,
                     cex_main = 15,
                     cex_axislab = 13,
                     cex_axistick = 13,
                     xmin = NULL,
                     xmax = NULL,
                     generatingnormal = NULL,
                     generatingamoroso = NULL,
                     generatingexgauss = NULL) {
  
  # Create main folder if it doesn't exist
  if (!dir.exists(newfilename)) {
    dir.create(newfilename)
  }
  
  # Create subfolder for raw images
  rawimg_folder <- file.path(newfilename, "rawimg")
  if (!dir.exists(rawimg_folder)) {
    dir.create(rawimg_folder)
  }
  
  # Calculate nr of total batches and size of partial batch
  num_full_batches <- length(dat) %/% batch_size  # Number of full 10-sized batches
  remaining <- length(dat) %% batch_size  # Remaining observations
  
  # Define max_y
  max_y <- max_y
  
  # Iteratively add 10 new data points each time
  for (i in 1:num_full_batches) {
    batch_data <- dat[1:(i * batch_size)]  # Take first i*10 observations
    glimpse(batch_data)
    
    tmp <- tempfile()
    svglite::svglite(tmp, width = svgwidth, height = svgheight)
    res <- estimate_methods(batch_data)
    
    if (!is.null(generatingnormal)) {
      pars <- generatingamoroso
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax, cex_main = cex_main,
                   cex_axislab = cex_axislab, cex_axistick = cex_axistick,
                   generatingnormal = c(pars[1],pars[2]))
    } else if (!is.null(generatingamoroso)) {
      pars <- generatingamoroso
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax, cex_main = cex_main,
                   cex_axislab = cex_axislab, cex_axistick = cex_axistick,
                   generatingamoroso = c(pars[1],pars[2],pars[3],pars[4]))
    } else if (!is.null(generatingexgauss)) {
      pars <- generatingexgauss
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax, cex_main = cex_main,
                   cex_axislab = cex_axislab, cex_axistick = cex_axistick,
                   generatingexgauss = c(pars[1],pars[2],pars[3]))
    } else {
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax, cex_main = cex_main,
                   cex_axislab = cex_axislab, cex_axistick = cex_axistick)
    }
    
    dev.off()
    
    
    # Construct dynamic filename
    filename <- file.path(rawimg_folder, sprintf("img%03d.png", i))
    # Save image
    rsvg_png(tmp, filename, width = pngwidth)
  }
  
  # Add remaining observations in the last batch
  if (remaining > 0) {
    batch_data <- dat  # Take all observations
    glimpse(batch_data)
    
    tmp <- tempfile()
    svglite::svglite(tmp, width = svgwidth, height = svgheight)
    res <- estimate_methods(batch_data)
    
    if (!is.null(generatingnormal)) {
      pars <- generatingamoroso
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax,
                   generatingnormal = c(pars[1],pars[2]))
    } else if (!is.null(generatingamoroso)) {
      pars <- generatingamoroso
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax,
                   generatingamoroso = c(pars[1],pars[2],pars[3],pars[4]))
    } else if (!is.null(generatingexgauss)) {
      pars <- generatingexgauss
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax,
                   generatingexgauss = c(pars[1],pars[2],pars[3]))
    } else {
      plot_methods(batch_data, res, ymax = max_y, yticks = c(0,max_y),
                   xmin = xmin, xmax = xmax)
    }
    
    dev.off()
    
    # Construct dynamic filename
    filename <- file.path(rawimg_folder, sprintf("img%03d.png", i))
    # Save image
    rsvg_png(tmp, filename, width = pngwidth)
    
  }
  
  print(paste("Plots saved in", newfilename, "folder."))
  
  # Create GIF from images
  imgs <- list.files(rawimg_folder, full.names = TRUE)
  img_list <- lapply(imgs, image_read)
  
  # Join and animate images
  img_joined <- image_join(img_list)
  img_animated <- image_animate(img_joined, fps = 2)
  
  # View animated image
  #img_animated
  
  # Save final GIF dynamically
  gif_filename <- file.path(paste0(newfilename,"/",newfilename,".gif"))
  image_write(image = img_animated, path = gif_filename)
  
  print(paste("GIF saved at:", gif_filename))
  
}

#make_gif()

