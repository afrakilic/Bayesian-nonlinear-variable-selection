#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######             MCMC MODEL SEARCH METHOD PERFORMANCE SIMULATION                  #####
#########################################################################################

source("R-codes/bayesian_selection.R")
source("R-codes/utils.R")

# Libraries 
library(mgcv)  # For generalized additive models
library(data.table)

#########################################################################################

# STATISTICAL BEHAVIOR ACROSS EFFECT SIZE

# Parameters
beta_values <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3)
sample_sizes <- c(100, 250)

results_beta <- list()
max_retries <- 10 # Set maximum retries for failed operations due to the convergence of gam. 

# Run Bayesian selection for each combination of beta and sample size
for (beta in beta_values) {
  for (n in sample_sizes) {
    key <- paste0("beta", beta, "_n_", n)
    retries <- 0
    
    repeat {
      retries <- retries + 1
      tryCatch(
        {
          results_beta[[key]] <- run_bayesian_selection(beta = beta, sample_size = n)
          break  # Exit loop on success
        },
        error = function(e) {
          message(sprintf("Error for beta = %.1f, n = %d (Retry %d): %s", beta, n, retries, e$message))
        }
      )
      if (retries >= max_retries) {
        message(sprintf("Failed for beta = %.1f, n = %d after %d retries.", beta, n, max_retries))
        break
      }
    }
  }
}


#Define trials 
trials <- 10

# Generalized function to compute statistics for given beta and column index
compute_stats <- function(sample_size, col_idx, stat_func = mean) {
  sapply(beta_values, function(beta) {
    key <- paste0("beta", beta, "_n_", sample_size)
    stat_func(results_beta[[key]][, col_idx])  # Apply statistic function
  })
}

# Correct selection ratio
cs_beta1 <- compute_stats(100, col_idx = 1, stat_func = function(x) sum(x == 1) / trials)
cs_beta2 <- compute_stats(250, col_idx = 1, stat_func = function(x) sum(x == 1) / trials)

# Posterior probabilities of the true model 
pp_beta1 <- compute_stats(100, col_idx = 2, stat_func = mean)
pp_beta2 <- compute_stats(250, col_idx = 2, stat_func = mean)


#PLOTS 
    
# Correct selection ratio
plot(beta_values, cs_beta1, type = "l", ylim = c(0, 1), 
     ylab = "Proportion of correct selection", xlab = "Sample Size")
lines(beta_values, cs_beta2, lty = 2)
legend("topright", legend = c(expression(paste(beta, " = 0.5")), expression(paste(beta, " = 2"))), lty = 1:2)

# Posterior probabilities of the true model
plot(beta_values, pp_beta1, type = "l", ylim = c(0, 1), 
     ylab = "Posterior Probability of the true model", xlab = "Sample Size")
lines(beta_values, pp_beta2, lty = 2)
legend("topright", legend = c(expression(paste(beta, " = 0.5")), expression(paste(beta, " = 2"))), lty = 1:2)

