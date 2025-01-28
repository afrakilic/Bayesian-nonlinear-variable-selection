

# Import CSV file
source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
source("R-codes/utils.R")

prelim.data <- read.csv("/Users/hakilic/Downloads/superconductivty+data (1)/train.csv", header = TRUE)
prelim.data$critical_temp <- prelim.data$critical_temp^(1/3)
y <- prelim.data$critical_temp
X <- prelim.data[, -length(data)]
# Combine input and output
X$Outcome <- y


# Data Split
set.seed(42)  
train_indices <- sample(1:nrow(X), size = 0.25 * nrow(X))

# training and testing
train_input <- X[train_indices, -ncol(X)]  # Input variables for training
train_output <- X[train_indices, ncol(X)]  # Output variable for training
test_input <- X[-train_indices, -ncol(X)]  # Input variables for testing
test_output <- X[-train_indices, ncol(X)]  # Output variable for testing

# Training
superconductivity_results <- bayesian_selection(X = train_input, y=train_output, knots=6, iteration = 2000)

selected_model <- superconductivity_results$`selected model`

results <- evaluate_model(
  selected_model =selected_model,
  train_input = train_input,
  train_output = train_output,
  test_input = test_input,
  test_output = test_output,
  knots = 10
)

