source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoThesis/refs/",
              "heads/main/estimate_amoroso.R"))


# What happens to Amoroso if I only add the min to the training data?
# (instead of min AND max)

set.seed(43)
dat <- rnorm(100)
min_dat <- min(dat)
max_dat <- max(dat)

train <- dat[1:50]
test <- dat[51:100]

minintrain <- min_dat %in% train
maxintrain <- max_dat %in% train

if(minintrain) {print("min in train")}
if(maxintrain) {print("max in train")}


#################################################
### Test what happens if only min is in train ###
#################################################

# Range of train and test data
print(range(train))
print(range(test)) 

# If min is not in train, remove it from test and add it to train
if(!minintrain) {
  test <- test[test != min_dat]  # Remove min from test
  train <- c(train, min_dat)  # Add min to train
  minintrain <- TRUE
}

# If max is in train, remove it from train and add it to test
if(maxintrain) {
  train <- train[train != max_dat]  # Remove max from train
  test <- c(test, max_dat)  # Add max to test
  maxintrain <- TRUE
}

# Updated range of train and test data
print(range(train))
print(range(test))


# Estimate Amorosos
res <- estimate_amoroso(train, plot=2, print_results = TRUE)

pars <- res$all_models %>%
  filter(method_ID == "WASS-CDF", space == "+") %>%
  select(a, l, c, mu) %>%
  as.numeric()

print(pars)

# Make x range way out of bounds to check if Amoroso makes zero predictions
xrange <- seq(-10,10,length.out=100)

# Make predictions from the extracted Amoroso fit
pred <- dgg4(x, pars[1], pars[2], pars[3], pars[4])

# Turn off scientific notation
options(scipen = 999)

# Create a two-row dataframe with rounded values
df <- data.frame(
  x = round(xrange, 3),
  pred = round(pred, 3)
)

# Print the dataframe
print(df)


######################################################
### Test what happens if min AND max are in train ###
######################################################

# Range of train and test data
print(range(train))
print(range(test)) 

# If min is not in train, remove it from test and add it to train
if(!minintrain) {
  test <- test[test != min_dat]  # Remove min from test
  train <- c(train, min_dat)  # Add min to train
  minintrain <- TRUE
}

# If max is in train, remove it from train and add it to test
if(!maxintrain) {
  test <- test[test != max_dat] # Remove max from test
  train <- c(train, max_dat)  # Add max to train
  maxintrain <- TRUE
}

# Updated range of train and test data
print(range(train))
print(range(test))

# Estimate Amorosos
res <- estimate_amoroso(train, plot=2, print_results = TRUE)

pars <- res$all_models %>%
  filter(method_ID == "WASS-CDF", space == "-") %>%
  select(a, l, c, mu) %>%
  as.numeric()

print(pars)

# Make x range way out of bounds to check if Amoroso makes zero predictions
xrange <- seq(-10,10,length.out=100)

# Make predictions from the extracted Amoroso fit
pred <- dgg4(x, pars[1], pars[2], pars[3], pars[4])

# Turn off scientific notation
options(scipen = 999)

# Create a two-row dataframe with rounded values
df <- data.frame(
  x = round(xrange, 3),
  pred = round(pred, 3)
)

# Print the dataframe
print(df)