
#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                         BOSTON HOUSING DATA                                 #####
#########################################################################################


source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
source("R-codes/utils.R") 
library(glmnet)

# Data preparation
data(boston, package = "spdep")

# X Matrix
X <- as.data.frame(cbind(
  "CRIM" = boston.c$CRIM,
  "ZN" = boston.c$ZN,
  "INDUS" = boston.c$INDUS,
  "CHAS" = boston.c$CHAS,
  "NOX" = boston.c$NOX,
  "RM" = boston.c$RM,
  "AGE" = boston.c$AGE,
  "DIS" = boston.c$DIS,
  "RAD" = boston.c$RAD,
  "TAX" = boston.c$TAX,
  "PTRATIO" = boston.c$PTRATIO,
  "B" = boston.c$B,
  "LSTAT" = boston.c$LSTAT
))

# Vector of responses
y <- boston.c$CMEDV
y <- y - mean(y) #centering y 

#Quadratic terms  COMMENT THIS PART WHEN FITTING BVS
for (i in 1:ncol(X)) {
  X[[paste0("V", i, "_quadratic")]] <- X[[i]]^2
}

# Interaction terms
variable_names <- colnames(X[, 1:13])
for (i in 1:(13 - 1)) {
  for (j in (i + 1):13) {
    X[[paste0(variable_names[i], "_x_", variable_names[j])]] <- X[[i]] * X[[j]]
  }
}

# Combine input and output
X$Outcome <- y

# Data Split
set.seed(42)  
train_indices <- sample(1:nrow(X), size = 0.75 * nrow(X))

# training and testing
train_input <- X[train_indices, -ncol(X)]  # Input variables for training
train_output <- X[train_indices, ncol(X)]  # Output variable for training
test_input <- X[-train_indices, -ncol(X)]  # Input variables for testing
test_output <- X[-train_indices, ncol(X)]  # Output variable for testing

# Training
boston_results <- bayesian_selection(X = train_input, y=train_output, knots=8, iteration = 2000)

selected_model <- boston_results$`selected model`

# Model Training with the Selected Model 
results <- evaluate_model(
  selected_model =selected_model,
  train_input = train_input,
  train_output = train_output,
  test_input = test_input,
  test_output = test_output,
  knots = 10
)
