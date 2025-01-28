

#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                                SML DATA                                      #####
#########################################################################################


source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
library(glmnet)

prelim.data <- read.table("/Users/hakilic/Downloads/NEW-DATA-2.T15.txt",header = TRUE,sep = " ")

prelim.data <- prelim.data[ , -c(1:2,19:21)]
colnum <- ncol(prelim.data)
prelim.data <- prelim.data[ , c(3:colnum,1)]
prelim.data$X24.Day_Of_Week <- round(prelim.data$X24.Day_Of_Week,digits = 0)
prelim.data$X24.Day_Of_Week <- factor(prelim.data$X24.Day_Of_Week)

dummy_variables <- as.data.frame(model.matrix(~ X24.Day_Of_Week  - 1, data = prelim.data))

prelim.data <- prelim.data[, !colnames(prelim.data) %in% "X24.Day_Of_Week"]
prelim.data <- cbind(prelim.data[, -ncol(prelim.data)], dummy_variables, prelim.data[, ncol(prelim.data), drop = FALSE])

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
SML_results <- bayesian_selection(X = train_input, y=train_output, knots=8, iteration = 2000)

selected_model <- SML_results$`selected model`

results <- evaluate_model(
  selected_model =selected_model,
  train_input = train_input,
  train_output = train_output,
  test_input = test_input,
  test_output = test_output,
  knots = 8
)





