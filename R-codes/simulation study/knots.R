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

# Run Bayesian selection for each knots value
for (knots in knots_values) {
  # Create a unique name for each knots value
  key <- paste0("k", knots)
  
  # Store results in the list
  results_knots[[key]] <- run_bayesian_selection(knots = knots)
}


######################################################################################################