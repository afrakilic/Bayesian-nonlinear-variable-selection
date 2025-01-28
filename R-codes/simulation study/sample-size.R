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

#STATISTICAL BEHAVIOR ACROSS SAMPLE SIZE 


# Parameters
sample_sizes <- c(50, 100, 250, 500)  # List of sample sizes
beta_values <- c(0.5, 2)  # List of beta values

results_n <- list()
max_retries <- 10  # Set maximum retries for failed operations due to the convergence of gam. 

# Loop over sample sizes and beta values
for (n in sample_sizes) {
  for (beta in beta_values) {
    # Create a unique key for the combination
    key <- paste0("n_", n, "_beta_", beta)
    retries <- 0
    success <- FALSE
    
    while (!success && retries < max_retries) {
      retries <- retries + 1
      tryCatch(
        {
          # Run Bayesian selection and store the results
          results_n[[key]] <- run_bayesian_selection(beta = beta, sample_size = n)
          success <- TRUE  # Mark as successful if no error occurs
        },
        error = function(e) {
          # Log the error and retry if necessary
          message(sprintf(
            "Error encountered for n = %d, beta = %.2f (Retry %d): %s", 
            n, beta, retries, conditionMessage(e)
          ))
        }
      )
    }
    
    if (!success) {
      # Log failure after exhausting retries
      message(sprintf(
        "Failed to process n = %d, beta = %.2f after %d retries.", 
        n, beta, max_retries
      ))
    }
  }
}


# Define trials
trials <- 100

# Generalized function to compute statistics for a given beta and column index
compute_stats <- function(beta, col_idx, stat_func = mean) {
  sapply(sample_sizes, function(n) {
    key <- paste0("n_", n, "_beta_", beta)
    stat_func(results_n[[key]][, col_idx])  # Apply statistic function (e.g., mean, sum)
  })
}

# Compute proportions of correct selection for beta = 0.5 and beta = 2
cs_n_0.5 <- compute_stats(0.5, col_idx = 1, stat_func = function(x) sum(x == 1) / trials)
cs_n_2 <- compute_stats(2, col_idx = 1, stat_func = function(x) sum(x == 1) / trials)

# Compute posterior probabilities of the true model for beta = 0.5 and beta = 2
pp_n_0.5 <- compute_stats(0.5, col_idx = 2, stat_func = mean)
pp_n_2 <- compute_stats(2, col_idx = 2, stat_func = mean)


#PLOTS 

# Plot proportion of correct selection
plot(sample_sizes, cs_n_0.5, type = "l", ylim = c(0, 1), 
     ylab = "Proportion of Correct Selection", xlab = "Sample Size")
lines(sample_sizes, cs_n_2, lty = 2)
legend("topright", legend = c(expression(paste(beta, " = 0.5")), expression(paste(beta, " = 2"))), lty = 1:2)

# Plot posterior probability of the true model
plot(sample_sizes, pp_n_0.5, type = "l", ylim = c(0, 1), 
     ylab = "Posterior Probability of the True Model", xlab = "Sample Size")
lines(sample_sizes, pp_n_2, lty = 2)
legend("topright", legend = c(expression(paste(beta, " = 0.5")), expression(paste(beta, " = 2"))), lty = 1:2)
