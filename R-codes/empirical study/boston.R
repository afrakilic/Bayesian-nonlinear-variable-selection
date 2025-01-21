
#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                         BOSTON HOUSING DATA                                 #####
#########################################################################################


source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
library(glmnet)
# Data preparation
data(boston, package = "spdep")

# X Matrix
X <- as.data.frame(cbind(
  "CRIM" = boston.c$CRIM,
  "ZN" = boston.c$ZN,
  "INDUS" = boston.c$INDUS,
  "CHAS" = boston.c$CHAS,
  "NOX" = boston.c$NOX,
  "RM" = boston.c$RM,
  "AGE" = boston.c$AGE,
  "DIS" = boston.c$DIS,
  "RAD" = boston.c$RAD,
  "TAX" = boston.c$TAX,
  "PTRATIO" = boston.c$PTRATIO,
  "B" = boston.c$B,
  "LSTAT" = boston.c$LSTAT
))

# Vector of responses
y <- boston.c$CMEDV
y <- y - mean(y)
#Quadratic terms
#for (i in 1:ncol(X)) {
 # X[[paste0("V", i, "_quadratic")]] <- X[[i]]^2
#}

# Interaction terms
variable_names <- colnames(X[, 1:13])
for (i in 1:(13 - 1)) {
  for (j in (i + 1):13) {
    X[[paste0(variable_names[i], "_x_", variable_names[j])]] <- X[[i]] * X[[j]]
  }
}

# Combine input and output
X$Outcome <- y

# Data Split
set.seed(42)  
train_indices <- sample(1:nrow(X), size = 0.75 * nrow(X))

# training and testing
train_input <- X[train_indices, -ncol(X)]  # Input variables for training
train_output <- X[train_indices, ncol(X)]  # Output variable for training
test_input <- X[-train_indices, -ncol(X)]  # Input variables for testing
test_output <- X[-train_indices, ncol(X)]  # Output variable for testing

# Training
boston_results <- bayesian_selection(X = train_input, y=train_output, knots=8, iteration = 2000)

selected_model <- boston_results$`selected model`

boston_91 <- boston_results

# Model Training with the Selected Model 

results <- selected_model_train(
  data = train_input,
  output = train_output,
  selected_model = selected_model,
  knots = 8,
  plot_title = "Training Predictions vs Real Output"
)

# Prediction 
pred_test <- predict.gam(results$model, test_input[, selected_model != 0])

results <- data.frame(
  Real_Output = test_output,
  Predicted_Output = pred_test
)

# Plot for the fitted values 
ggplot(results, aes(x = Real_Output, y = Predicted_Output)) +
  geom_point(color = "blue", alpha = 0.6) + # Scatter points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + # Line y=x
  labs(
    title = "Predictions vs Real Output",
    x = "Real Output (y)",
    y = "Predicted Output"
  ) +
  theme_minimal()

rmse_test <- sqrt(mean((test_output - pred_test)^2))

print(paste("RMSE:", rmse_test))



#LASSO 
x = as.matrix(train_input)
y = train_output
cv_lasso <- cv.glmnet(x, y, alpha = 1)

# Best lambda value
best_lambda <- cv_lasso$lambda.min
print(best_lambda)
final_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)

# Predict using the final model
predictions <- predict(final_model, s = best_lambda, newx = as.matrix(test_input))
rmse_test_lasso <- sqrt(mean((test_output - predictions)^2))
rmse_test_lasso





  
  
  
  
  
  
  
  
  

