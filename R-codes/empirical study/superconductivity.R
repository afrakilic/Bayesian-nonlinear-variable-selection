

# Import CSV file
data <- read.csv("/Users/hakilic/Downloads/superconductivty+data (1)/train.csv", header = TRUE)

y <- data$critical_temp
X <- data[, -length(data)]
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

knots = 8
selected_model <- boston_91$`selected model`
data = train_input
# for linear effects
linears <- c()
if(length(selected_model[selected_model == 1]) != 0){
  
  for(i in 1:ncol(data[which(selected_model == 1)])){
    linears <- c(linears, paste(c(colnames(data[which(selected_model == 1)][i])), collapse= ""))
  } 
}

#for nonlinear effects
non_linears <- c()
if(length(selected_model[selected_model == 2]) != 0) {
  for(i in 1:ncol(data[which(selected_model == 2)])){
    non_linears <- c(non_linears, paste(c('s(', colnames(data[which(selected_model == 2)][i]), ',k=',knots,')'), collapse= ""))
  } 
}


vars <- c(linears, non_linears)
y <- train_output 

model = gam(as.formula(paste('y', '~ 1 + ', paste(vars, collapse =  "+"))), data = data)
preds = predict(model, data)
# Combine predictions and actual values into a dataframe for easy plotting
results <- data.frame(
  Real_Output = y,
  Predicted_Output = preds
)

# Create a scatter plot to compare predictions vs real output
ggplot(results, aes(x = Real_Output, y = Predicted_Output)) +
  geom_point(color = "blue", alpha = 0.6) + # Scatter points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + # Line y=x
  labs(
    title = "Predictions vs Real Output",
    x = "Real Output (y)",
    y = "Predicted Output"
  ) +
  theme_minimal()

# Calculate RMSE
rmse_train <- sqrt(mean((y - preds)^2))
print(paste("RMSE:", rmse_train))

# Prediction 
# Prediction 
pred_test <- predict.gam(model, test_input[, selected_model != 0])

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


