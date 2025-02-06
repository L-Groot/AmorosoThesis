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

# Glimpse
glimpse(acidity)
glimpse(chondrite)
glimpse(enzyme)
glimpse(galaxy)
glimpse(geyser)
glimpse(stamps)

# Fit the 5 methods to each of the 6 datasets
estimate_methods(acidity, amoinaplus = F, minimal = T)
estimate_methods(chondrite, amoinaplus = F, minimal = T)
estimate_methods(enzyme, amoinaplus = F, minimal = T)
estimate_methods(galaxy, amoinaplus = F, minimal = T, yticks = c(0.00000,0.00025))
estimate_methods(geyser, amoinaplus = F, minimal = T)
estimate_methods(stamps, amoinaplus = F, minimal = T)

glimpse(chondrite)
