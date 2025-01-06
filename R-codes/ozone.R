#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                                OZONE DATA                                   #####
#########################################################################################

source("bayesian_selection.R")  # Ensure the script is in the working directory

#OZONE DATA 
library(mlbench)
data(Ozone)

# Prepare the data
Ozone <- Ozone[, -c(1, 2, 3, 9)]  # Remove specified columns
y_ozone <- Ozone$V4               # Extract outcome variable
Ozone <- Ozone[, -1]              # Remove the outcome variable

# Rename columns sequentially as Variable1, Variable2, ...
colnames(Ozone) <- paste0("Variable", 1:ncol(Ozone))

# Square terms
for (i in 1:ncol(Ozone)) {
  Ozone[[paste0("V", i, "_square")]] <- Ozone[[i]]^2
}

# Interaction terms
variable_names <- colnames(Ozone[, 1:8])
for (i in 1:(8 - 1)) {
  for (j in (i + 1):8) {
    Ozone[[paste0(variable_names[i], "_x_", variable_names[j])]] <- Ozone[[i]] * Ozone[[j]]
  }
}

# Combine input and output
Ozone$Outcome <- y_ozone

# Split into train (75%) and test (25%)
set.seed(42)  # For reproducibility
train_indices <- sample(1:nrow(Ozone), size = 0.75 * nrow(Ozone))

# Separate training and testing data
train_input <- Ozone[train_indices, -ncol(Ozone)]  # Input variables for training
train_output <- Ozone[train_indices, ncol(Ozone)] # Output variable for training
test_input <- Ozone[-train_indices, -ncol(Ozone)]  # Input variables for testing
test_output <- Ozone[-train_indices, ncol(Ozone)] # Output variable for testing


# Training
ozone_results <- bayesian_selection(data_original = train_input, knots = 8, y= train_output, iteration = 2000)

ozone_results$`selected model`
ozone_results$`posterior probability_selected`

# Prediction

#########################################################################################