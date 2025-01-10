#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                             DIABETIES DATA                                   #####
#########################################################################################

source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory

data(diabetesI, package = "spikeslab")

y <- diabetesI$Y #the outcome
Diabeties <- diabetesI[, -1] # removing the outcome 

# Combine input and output
Diabeties$Outcome <- y

# Split into train (75%) and test (25%)
set.seed(42)  # For reproducibility
train_indices <- sample(1:nrow(Diabeties), size = 0.75 * nrow(Diabeties))

# Separate training and testing data
train_input <- Diabeties[train_indices, -ncol(Diabeties)]  # Input variables for training
train_output <- Diabeties[train_indices, ncol(Diabeties)] # Output variable for training
test_input <- Diabeties[-train_indices, -ncol(Diabeties)]  # Input variables for testing
test_output <- Diabeties[-train_indices, ncol(Diabeties)] # Output variable for testing

# Training
diabetes_results <- bayesian_selection(X = train_input, knots = 8, y= train_output, iteration = 2000)


diabetes_results$`selected model`
diabetes_results$`posterior probability_selected`

# Model Training with the Selected Model 

knots = 8
selected_model <- diabetes_results$`selected model`
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
pred_train = predict(model, data)
# Combine predictions and actual values into a dataframe for easy plotting
results <- data.frame(
  Real_Output = y,
  Predicted_Output = pred_train
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

# RMSE
rmse_train <- sqrt(mean((y - pred_train)^2))
print(paste("RMSE:", rmse_train))

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



