#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                                OZONE DATA                                   #####
#########################################################################################

source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
source("R-codes/bayesian_selection_large.R")
#OZONE DATA 
library(mlbench)
data(Ozone)

# Prepare the data
Ozone <- Ozone[, -c(1, 2, 3, 9)]
Ozone <- apply(Ozone, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})
Ozone <- as.data.frame(Ozone)
# Remove specified columns
y <- Ozone$V4               # Extract outcome variable
Ozone <- Ozone[, -1]              # Remove the outcome variable

# Rename columns sequentially as Variable1, Variable2, ...
colnames(Ozone) <- paste0("Variable", 1:ncol(Ozone))

 #Square terms
for (i in 1:ncol(Ozone)) {
  Ozone[[paste0("V", i, "_square")]] <- Ozone[[i]]^2
}

# Interaction terms
variable_names <- colnames(Ozone[, 1:8])
for (i in 1:(8 - 1)) {
  for (j in (i + 1):8) {
    Ozone[[paste0(variable_names[i], "_x_", variable_names[j])]] <- Ozone[[i]] * Ozone[[j]]
  }
}

# Combine input and output
Ozone$Outcome <- y

# Split into train (75%) and test (25%)
set.seed(42)  # For reproducibility
train_indices <- sample(1:nrow(Ozone), size = 0.75 * nrow(Ozone))

# Separate training and testing data
train_input <- Ozone[train_indices, -ncol(Ozone)]  # Input variables for training
train_output <- Ozone[train_indices, ncol(Ozone)] # Output variable for training
test_input <- Ozone[-train_indices, -ncol(Ozone)]  # Input variables for testing
test_output <- Ozone[-train_indices, ncol(Ozone)] # Output variable for testing


# Training
ozone_results <- bayesian_selection_large(X = train_input, knots = 8, y= train_output, iteration = 2000)

ozone_results$`selected model`
ozone_results$`posterior probability_selected`

# Model Training with the Selected Model 

knots = 8
selected_model <- ozone_results$`selected model`
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
