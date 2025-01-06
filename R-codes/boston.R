
#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                         BOSTON HOUSING DATA                                 #####
#########################################################################################


source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory

# Data preparation
library(spdep)
data(boston)

# X Matrix
X_boston <- as.data.frame(cbind(
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
y_boston <- boston.c$CMEDV

# Quadratic terms
for (i in 1:ncol(X_boston)) {
  X_boston[[paste0("V", i, "_quadratic")]] <- X_boston[[i]]^2
}

# Interaction terms
variable_names <- colnames(X_boston[, 1:13])
for (i in 1:(13 - 1)) {
  for (j in (i + 1):13) {
    X_boston[[paste0(variable_names[i], "_x_", variable_names[j])]] <- X_boston[[i]] * X_boston[[j]]
  }
}

# Combine input and output
X_boston$Outcome <- y_boston

# Data Split
set.seed(42)  
train_indices <- sample(1:nrow(X_boston), size = 0.75 * nrow(X_boston))

# training and testing
train_input <- X_boston[train_indices, -ncol(X_boston)]  # Input variables for training
train_output <- X_boston[train_indices, ncol(X_boston)]  # Output variable for training
test_input <- X_boston[-train_indices, -ncol(X_boston)]  # Input variables for testing
test_output <- X_boston[-train_indices, ncol(X_boston)]  # Output variable for testing

# Training
boston_results <- bayesian_selection(data_original = train_input, y=train_output, knots=6, iteration = 2000)

boston_results$`selected model`

bostonn <- boston_results

# Prediction
#########################################################################################


