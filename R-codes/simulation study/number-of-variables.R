#### Title: Bayesian Variable Selection for Linear and Nonlinear Models
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

# Function to run Bayesian selection and collect results for a given n_var value
run_bayesian_for_n_var <- function(n_var, trials = 100) {
  results <- matrix(NA, nrow = trials, ncol = 3)  # Pre-allocate matrix for results
  times <- numeric(trials)  # To store CPU times for each trial
  selected_model <- matrix(NA, nrow = trials, ncol = n_var)  # Pre-allocate matrix
  
  for (i in seq_len(trials)) {
    # Start the timer
    start <- proc.time()
    
    # Generate data and apply Bayesian selection
    data <- generate_true_model(n_var = n_var)
    model <- bayesian_selection(X = data$data, y = data$response)
    
    # Compute posterior probabilities
    posterior <- calculate_posterior_probability(model$`gamma draws`, data$true_model)
    results[i, ] <- c(
      posterior$`Is selected model the true model`,
      posterior$posterior_prob_true_model,
      posterior$posterior_prob_selected_model
    )
    selected_model[i, ] <- model$`selected model`
    
    # Record elapsed time
    times[i] <- proc.time()["elapsed"] - start["elapsed"]
  }
  
  # Return results, models, and timings
  list(results = results, gammas = selected_model, times = times)
}

# Define n_vars and run simulations
n_vars <- c(6, 10, 15, 20, 30)
results_list <- lapply(n_vars, run_bayesian_for_n_var)

# Extract results, CPU times, and gamma draws
results <- lapply(results_list, `[[`, "results")
cpu_times <- sapply(results_list, function(x) mean(x$times))
gamma_draws <- lapply(results_list, `[[`, "gammas")

# Compute summary table (J)
J <- do.call(cbind, lapply(results, function(res) {
  correct_selection <- mean(res[, 1])  # Proportion of correct selections
  posterior_prob <- mean(res[, 2])     # Mean posterior probability of the true model
  c(correct_selection, posterior_prob)
}))

# Add CPU time ratios to J
cpu_ratio <- cpu_times / cpu_times[2]
J <- rbind(J, cpu_ratio)
rownames(J) <- c("Correct Selection Proportion", "Posterior Pr. of the True Model", "CPU")
colnames(J) <- paste0("J=", n_vars)

# Display results
print(J)

#########################################################################################

# MULTIPLICITY CHECK

# Function to compute row-wise counts for 0, 1, and 2
count_values <- function(mat) {
  t(apply(mat, 1, function(row) tabulate(factor(row, levels = 0:2), nbins = 3)))
}

# Compute multiplicity
multiplicity <- do.call(rbind, lapply(gamma_draws, function(gammas) {
  colMeans(count_values(gammas))
}))
rownames(multiplicity) <- paste0("J=", n_vars)

# Display multiplicity results
print(multiplicity)

# Plot results
plot(n_vars, multiplicity[, 3], type = "o", col = "darkblue", ylim = c(0, max(multiplicity)),
     main = "Number of Identified Effect Types Across J", xlab = "J", ylab = "Number of Identified Effect Type")
lines(n_vars, multiplicity[, 2], type = "o", col = "black")
lines(n_vars, multiplicity[, 1], type = "s", col = "black", lty = 2)
points(n_vars, multiplicity[, 1], pch = 18)
legend("topright", legend = c(expression(paste(gamma, "=1")), expression(paste(gamma, "=2")), expression(paste(gamma, "=0"))),
       col = c("darkblue", "black", "black"), lty = c(1, 1, 2), pch = c(NA, NA, 18), cex = 0.8)

