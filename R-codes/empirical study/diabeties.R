#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                             DIABETIES DATA                                   #####
#########################################################################################

source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
source("R-codes/utils.R")
library(glmnet)

data(diabetesI, package = "spikeslab")

y <- diabetesI$Y #the outcome
y <- y - mean(y)
X <- diabetesI[, -1] # removing the outcome 
X <- X[, -c(11:19)] #removing the quadratics for gam (COMMENT THIS PART WHEN FITTING LASSO)
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
diabetes_results <- bayesian_selection(X = train_input, knots = 8, y= train_output, iteration = 2000)
selected_model = diabetes_results$`selected model`

results <- evaluate_model(
  selected_model =selected_model,
  train_input = train_input,
  train_output = train_output,
  test_input = test_input,
  test_output = test_output,
  knots = 8
)
