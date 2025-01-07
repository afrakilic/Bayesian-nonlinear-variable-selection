#### Title: Bayesian Variable Selection for Linear and Nonlinear Models
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                            HELPER FUNCTIONS                                  #####
#########################################################################################

# Load necessary libraries
library(mgcv)          # For generalized additive models
library(data.table)    # Efficient data manipulation
library(MCMCprecision) # Precision for Bayesian MCMC
library(lubridate)     # Date-time operations
library(progress)      # Progress bar for iterative processes

#########################################################################################
######                  FUNCTION: generate_true_model                               #####
#########################################################################################

# Function to simulate the true model for Bayesian variable selection
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

# Function to calculate posterior probabilities of the true and selected models
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


run_bayesian_selection <- function(beta = 0.8, sample_size = 100, n_var = 10, trials = 10, knots = 4) {
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


