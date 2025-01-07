######################################################################################################
#EFFECT SIZES (BETA)
######################################################################################################

# Define parameters
beta_values <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3)
sample_sizes <- c(100, 250)
results_beta <- list()
max_retries <- 10

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
######################################################################################################
# Define sample sizes and trials
trials <- 10

# Generalized function to compute statistics for given beta and column index
compute_stats <- function(beta, col_idx, stat_func = mean) {
  sapply(sample_sizes, function(n) {
    key <- paste0("beta_", beta, "_n_", n)
    stat_func(results_n[[key]][, col_idx])  # Apply statistic function (e.g., mean, sum)
  })
}

# Compute proportions of correct selection for beta = 0.5 and beta = 2
cs_beta_100 <- compute_stats(0.5, col_idx = 1, stat_func = function(x) sum(x == 1) / trials)
cs_n_2 <- compute_stats(2, col_idx = 1, stat_func = function(x) sum(x == 1) / trials)

# Compute posterior probabilities of the true model for beta = 0.5 and beta = 2
pp_n_0.5 <- compute_stats(0.5, col_idx = 2, stat_func = mean)
pp_n_2 <- compute_stats(2, col_idx = 2, stat_func = mean)

# Function to plot results
plot_results <- function(x, y1, y2, y_label, legend_labels, ylim = c(0, 1), legend_pos = "topright") {
  plot(x, y1, type = "l", ylim = ylim, ylab = y_label, xlab = "Sample Size")
  lines(x, y2, type = "l", lty = 2)
  legend(legend_pos, legend = legend_labels, lty = 1:2, cex = 1)
}

# Plot proportion of correct selection
plot_results(
  sample_sizes, cs_n_0.5, cs_n_2, 
  y_label = "Proportion of correct selection", 
  legend_labels = c(expression(paste(beta, " = 0.5")), expression(paste(beta, " = 2")))
)

# Plot posterior probability of the true model
plot_results(
  sample_sizes, pp_n_0.5, pp_n_2, 
  y_label = "Posterior Probability of the true model", 
  legend_labels = c(expression(paste(beta, " = 0.5")), expression(paste(beta, " = 2")))
)
