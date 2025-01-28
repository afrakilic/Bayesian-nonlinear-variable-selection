
#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                               COLLEGE DATA                                  #####
#########################################################################################


source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
library(glmnet)

data("College", package= "ISLR")
prelim.data=College

col.num <- ncol(prelim.data)
prelim.data=prelim.data[ , c(1,3:(col.num),2)]
# removed enroll and accept due to causual issues
prelim.data=subset(prelim.data, select=-c(Enroll,Accept))
prelim.data$Apps=log(prelim.data$Apps)
prelim.data$F.Undergrad=log(prelim.data$F.Undergrad)
prelim.data$P.Undergrad=log(prelim.data$P.Undergrad)
prelim.data[, 1] <- ifelse(prelim.data[, 1] == "Yes", 1, 0)
# centering the Y. Assuming Y is always the last column of prelim.data
y.val.0 <- mean(prelim.data[ , ncol(prelim.data)])
prelim.data[ , ncol(prelim.data)]=prelim.data[ , ncol(prelim.data)]-mean(prelim.data[ , ncol(prelim.data)])

X <- prelim.data

# Data Split
set.seed(42)  
train_indices <- sample(1:nrow(X), size = 0.75 * nrow(X))

# training and testing
train_input <- X[train_indices, -ncol(X)]  # Input variables for training
train_output <- X[train_indices, ncol(X)]  # Output variable for training
test_input <- X[-train_indices, -ncol(X)]  # Input variables for testing
test_output <- X[-train_indices, ncol(X)]  # Output variable for testing

# Training
college_results <- bayesian_selection(X = train_input, y=train_output, knots=10, iteration = 2000)

selected_model <- college_results$`selected model`

results <- evaluate_model(
  selected_model =selected_model,
  train_input = train_input,
  train_output = train_output,
  test_input = test_input,
  test_output = test_output,
  knots = 10
)
