#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######             MCMC MODEL SEARCH METHOD PERFORMANCE SIMULATION                  #####
#########################################################################################

source("R-codes/bayesian_selection.R")
source("R-codes/utils.R")

# libraries 
library(mgcv)  # For generalized additive models
library(data.table)

#########################################################################################

# STATISTICAL BEHAVIOR ACROSS NUMBER OF BASIS FUNCTION USED IN MODEL TRAINING

knots_values <- c(4, 6, 8, 10, 12)

# Initialize an empty list to store results
results_knots <- list()

# Initialize a vector to store CPU time for each knots value
CPU <- numeric(length(knots_values))

# Run Bayesian selection for each knots value and record the computation time
for (i in seq_along(knots_values)) {
  knots <- knots_values[i]
  
  # Measure computation time for each knots value
  time_taken <- system.time({
    # Store results in the list
    key <- paste0("k", knots)
    results_knots[[key]] <- run_bayesian_selection(knots = knots)
  })
  
  # Record elapsed time in seconds
  CPU[i] <- time_taken["elapsed"]
}

# Statistics for given column index and knots values
compute_stats_knots <- function(col_idx, stat_func = mean) {
  sapply(knots_values, function(knots) {
    key <- paste0("k", knots)  # Create the key based on knots value
    stat_func(results_knots[[key]][, col_idx])  # Apply the statistic function to the specified column
  })
}

trials = 100
#Correct selection ratio
knots_correct_select <- compute_stats_knots(col_idx = 1, stat_func = function(x) sum(x == 1) / trials)
#Posterior probabilities of the true model 
knots_post_prob <- compute_stats_knots(col_idx = 2, stat_func = mean)

knots_correct_select
knots_post_prob


######################################################################################################