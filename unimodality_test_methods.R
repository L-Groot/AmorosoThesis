# Unimodality methods

library(tidyverse)
library(multimode)
library(BayesMultiMode)

# Load data 
affect <- read.csv("data_affect.csv")

# Take a glimpse at data
glimpse(affect)

# Create 8 subsets
emo1 <- affect$emo1_m
emo2 <- affect$emo2_m
emo3 <- affect$emo3_m
emo4 <- affect$emo4_m
emo5 <- affect$emo5_m
emo6 <- affect$emo6_m
emo7 <- affect$emo7_m
emo8 <- affect$emo8_m

# Split plotting window
par(mfrow = c(3,3))

# Plot histograms
hist(emo1, main = "emo1")
hist(emo2, main = "emo2")
hist(emo3, main = "emo3")
hist(emo4, main = "emo4")
hist(emo5, main = "emo5")
hist(emo6, main = "emo6")
hist(emo7, main = "emo7")
hist(emo8, main = "emo8")


# Bimodality coefficient
# Hartigan's d
# GMM
# Haslbecks method
# Silverman
# Excess-mass

# d1b1g1 <- affect %>% 
#   filter(dayno == 1) %>%
#   select(emo3_m) %>%
#   pull()


BayesMultiMode::

