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

# CONVERGENCE OF GIBSS SAMPLER CHAIN #

compute_convergence <- function(n_var) {
  data <- generate_true_model(n_var = n_var)
  chain <- bayesian_selection(X = data$data, y = data$response, iteration = 10000)
  
  chain_pp <- numeric(100)  # Matrix to store posterior probabilities
  for (i in 1:100) {
    gamma_draws <- chain$`gamma draws`
    posterior <- data.table(gamma_draws[1:(i * 100), ])
    frequency <- as.matrix(posterior[, list(posterior = .N), by = names(posterior)][order(posterior, decreasing = TRUE)])
    chain_pp[i] <- frequency[1, (n_var + 1)] / (i * 100)  # Posterior probability calculation
  }
  
  plot(100 * (1:100), chain_pp, type = "l", ylim = c(0, 1), 
       ylab = "Posterior probability of the selected model", xlab = "T", main = paste("J =", n_var))
}

# Plot for J=6 and J=20
compute_convergence(6)
compute_convergence(20)
