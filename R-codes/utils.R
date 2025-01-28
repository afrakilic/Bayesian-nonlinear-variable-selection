#### Title: Bayesian Variable Selection for Linear and Nonlinear Models
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                            HELPER FUNCTIONS                                  #####
#########################################################################################

# Libraries
library(mgcv)# For generalized additive models
library(data.table)   
library(ggplot2)
library(glmnet)
source("R-codes/bayesian_selection.R")  

# This script implements a simulation framework for Bayesian variable selection in linear 
# and nonlinear models. The goal is to identify which predictor variables have linear, 
# nonlinear, or no effect on the response variable by using Bayesian methods. The script 
# contains helper functions for generating synthetic data, calculating posterior probabilities, 
# and running simulations.


#########################################################################################
######                  FUNCTION: generate_true_model                               #####
#########################################################################################

# Function to simulate the true model for Bayesian variable selection
# This function generates synthetic data with a mix of linear, nonlinear, and zero-effect 
# predictors. It also creates the true model structure, which is later used for comparison 
# during Bayesian inference.

generate_true_model <- function(beta = 0.8, 
                                n = 100,    # Sample size
                                sigma2 = 0.1, # Variance of the error term
                                n_var = 10   # Number of predictor variables
) {
  # Data Generation
  # Create a matrix of independent predictor variables with random normal values
  variables <- matrix(NA, nrow = n, ncol = n_var)
  colnames(variables) <- paste0("x", 1:n_var) # Assign variable names x1, x2, ..., xn
  for (i in 1:n_var) {
    variables[, i] <- rnorm(n) # 
  }
  
  # Generate the outcome variable with random noise
  error <- rnorm(n, sd = sqrt(sigma2)) 
  y <- error #the response variable with random noise
  
  # true relationships in the model
  nonlinear_true <- c("x1", "x2") # Variables with nonlinear effects
  linear_true <- c("x3", "x4")   # Variables with linear effects
  zero_true <- setdiff(colnames(variables), c(nonlinear_true, linear_true)) # No effect variables
  
  # nonlinear effects 
  nonlinear_indices <- match(nonlinear_true, colnames(variables))
  for (i in seq_along(nonlinear_indices)) {
    y <- y + beta * exp(variables[, nonlinear_indices[i]])
  }
  
  # zero effects (illustrative, does not modify y)
  zero_indices <- match(zero_true, colnames(variables))
  for (i in seq_along(zero_indices)) {
    y <- y + 0 * variables[, zero_indices[i]] # No effect, included for clarity
  }
  
  # linear effects
  linear_indices <- match(linear_true, colnames(variables))
  for (i in seq_along(linear_indices)) {
    y <- y + beta * variables[, linear_indices[i]]
  }
  
  # Center the response variable to remove intercept
  y <- y - mean(y)
  
  # true model sequence
  true_gamma <- rep(NA, n_var) # Initialize gamma vector
  true_gamma[nonlinear_indices] <- 2 # Mark nonlinear effects as 2
  true_gamma[zero_indices] <- 0      # Mark zero effects as 0
  true_gamma[linear_indices] <- 1    # Mark linear effects as 1
  
  
  data_original <- as.data.frame(variables)
  
  result <- list(
    "data" = data_original,
    "response" = y,
    "true_model" = true_gamma
  )
  
  return(result)
}

#########################################################################################
######            FUNCTION: calculate_posterior_probability                        #####
#########################################################################################

# Function to calculate posterior probabilities of the true model and the most frequent model
# This function evaluates how well the Bayesian variable selection approach identifies the 
# true model structure by comparing posterior probabilities.

calculate_posterior_probability <- function(gamma_draws, true_gamma) {
  # Dimensions of the gamma draws matrix
  iteration <- dim(gamma_draws)[1] # Number of iterations
  n_var <- dim(gamma_draws)[2]    # Number of variables
  
  # Convert gamma draws to a data table for manipulation
  posterior <- data.table(gamma_draws)
  
  # Calculate frequency of each gamma sequence
  frequency <- as.matrix(
    posterior[, list(posterior = .N), by = names(posterior)][order(posterior, decreasing = TRUE)]
  )
  
  selected_model = as.vector(frequency[1,1:n_var])
  
  # Check if the true gamma sequence appears among the draws
  t <- apply(frequency[, -(n_var + 1)], 1, function(x) return(all(x == true_gamma)))
  
  # Calculate posterior probability of the true model
  pp_t <- if (any(t) == FALSE) 0 else as.numeric(frequency[which(t), (n_var + 1)]) / iteration
  
  # Calculate posterior probability of the selected (most frequent) model
  pp_s <- as.numeric(frequency[1, (n_var + 1)]) / iteration
  
  # Return posterior probabilities as a list
  result <- list(
    "posterior_prob_true_model" = pp_t,
    "posterior_prob_selected_model" = pp_s,
    "Is selected model the true model" = identical(true_gamma, selected_model)
  )
  
  return(result)
}


#########################################################################################
######            FUNCTION: running simulation for different models                 #####
#########################################################################################

# Function to run multiple trials of Bayesian variable selection
# This function repeatedly generates data, applies Bayesian selection, and evaluates
# the posterior probabilities of the true and selected models.

run_bayesian_selection <- function(beta = 0.8, sample_size = 100, n_var = 10, trials = 100, knots = 4) {
  results <- matrix(NA, nrow = trials, ncol = 3)  # Pre-allocate matrix for results
  
  for (i in 1:trials) {
    
    data <- generate_true_model(beta = beta, n = sample_size, n_var = n_var) 
    model <- bayesian_selection(X = data$data, y = data$response, knots = knots)
    a = calculate_posterior_probability(model$`gamma draws`, data$true_model)
    # Store results for each trial
    results[i, ] <- c(a$`Is selected model the true model`, a$posterior_prob_true_model, a$posterior_prob_selected_model)
  }
  
  return(results)
}

#########################################################################################
######            FUNCTION: empirical study BVS vs. LASSO                          #####
#########################################################################################

# This function evaluates and compares the performance of two models: 
# a user-selected model (which could include both linear and nonlinear components) 
# and a LASSO regression model. It calculates the inclusion rates of linear and 
# nonlinear variables for the selected model, computes the Root Mean Squared Error (RMSE) 
# for both the selected model and the LASSO model, and outputs a table comparing 
# these metrics. The function takes training and testing data, as well as the number of 
# knots to be used in the selected model, as inputs.

evaluate_model <- function(selected_model, train_input, train_output, test_input, test_output, knots = 3) {
  # Inclusion Rates
  nonlinear_rate <- sum(selected_model == 2) / ncol(train_input)
  linear_rate <- sum(selected_model == 1) / ncol(train_input)
  
  # Model Training with the Selected Model
  results <- selected_model_train(
    data = train_input,
    output = train_output,
    selected_model = selected_model,
    knots = knots,
    plot_title = "Training Predictions vs Real Output"
  )
  
  # Prediction
  pred_test <- predict.gam(results$model, test_input[, selected_model != 0])
  
  # RMSE for the selected model
  rmse_test <- sqrt(mean((test_output - pred_test)^2))
  
  # LASSO
  x <- as.matrix(train_input)
  y <- train_output
  cv_lasso <- cv.glmnet(x, y, alpha = 1)
  
  # Best lambda value
  best_lambda <- cv_lasso$lambda.min
  
  # Final LASSO model
  final_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  
  # Number of variables included
  total_included <- sum(coef(final_model) != 0) - 1
  inclusion_rate_lasso <- total_included / ncol(train_input)
  
  # Predict using the LASSO model
  predictions <- predict(final_model, s = best_lambda, newx = as.matrix(test_input))
  rmse_test_lasso <- sqrt(mean((test_output - predictions)^2))
  
  # Create a results table
  results_table <- data.frame(
    Metric = c(
      "Nonlinear Variable Inclusion Rate",
      "Linear Variable Inclusion Rate",
      "Selected Model RMSE",
      "LASSO Inclusion Rate",
      "LASSO RMSE"
    ),
    Value = c(
      nonlinear_rate,
      linear_rate,
      rmse_test,
      inclusion_rate_lasso,
      rmse_test_lasso
    )
  )
  
  # Print the results table
  print(results_table)
  
  # Return the results table
  return(results_table)
}
