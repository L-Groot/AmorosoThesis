# (1) Qualititative goodness of fit

# Use the 6 empirical datasets from Almeijeiras-Alonso et al. (2021) which are
# commonly used to evaluate the behavior of methods that aim to determine the
# number of modes

# Libraries
library(multimode)

# Load functions
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_methods.R"))

# Load the six datasets
data("acidity", package = "multimode")
data("chondrite", package = "multimode")
data("enzyme", package = "multimode")
data("galaxy", package = "multimode")
data("geyser", package = "multimode")
data("stamps", package = "multimode")

# Fit the 5 methods to each of the 6 datasets
estimate_methods(acidity, amoinaplus = FALSE)
estimate_methods(chondrite, breaks = 10, amoinaplus = FALSE)
estimate_methods(enzyme, amoinaplus = FALSE)
