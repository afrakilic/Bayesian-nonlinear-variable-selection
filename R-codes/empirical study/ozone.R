#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                                OZONE DATA                                   #####
#########################################################################################

source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
library(glmnet)

#OZONE DATA 
data(Ozone, package = "mlbench")

# Prepare the data
X <- Ozone[, -c(1, 2, 3, 9)]
X <- apply(X, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})
X <- as.data.frame(X)
y <- X$V4               # Extract outcome variable
y <- y - mean(y)
X <- X[, -1]              # Remove the outcome variable

colnames(X) <- paste0("Variable", 1:ncol(X))

 #Square terms COMMENT THIS PART WHEN FITTING BVS
for (i in 1:ncol(X)) {
  X[[paste0("V", i, "_square")]] <- X[[i]]^2
}

# Interaction terms
variable_names <- colnames(X[, 1:8])
for (i in 1:(8 - 1)) {
  for (j in (i + 1):8) {
    X[[paste0(variable_names[i], "_x_", variable_names[j])]] <- X[[i]] * X[[j]]
  }
}

# Combine input and output
X$Outcome <- y

# Split into train (75%) and test (25%)
set.seed(42)  # For reproducibility
train_indices <- sample(1:nrow(X), size = 0.75 * nrow(X))

# Separate training and testing data
train_input <- X[train_indices, -ncol(X)]  # Input variables for training
train_output <- X[train_indices, ncol(X)] # Output variable for training
test_input <- X[-train_indices, -ncol(X)]  # Input variables for testing
test_output <- X[-train_indices, ncol(X)] # Output variable for testing


# Training
ozone_results <- bayesian_selection(X = train_input, knots = 8, y= train_output, iteration = 2000)

selected_model <- ozone_results$`selected model`

results <- evaluate_model(
  selected_model =selected_model,
  train_input = train_input,
  train_output = train_output,
  test_input = test_input,
  test_output = test_output,
  knots = 8
)
