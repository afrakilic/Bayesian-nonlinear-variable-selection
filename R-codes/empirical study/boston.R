
#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                         BOSTON HOUSING DATA                                 #####
#########################################################################################


source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory

# Data preparation
library(spdep)
data(boston)

# X Matrix
X_boston <- as.data.frame(cbind(
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
y_boston <- boston.c$CMEDV

# Quadratic terms
#for (i in 1:ncol(X_boston)) {
#  X_boston[[paste0("V", i, "_quadratic")]] <- X_boston[[i]]^2
#}

# Interaction terms
variable_names <- colnames(X_boston[, 1:13])
for (i in 1:(13 - 1)) {
  for (j in (i + 1):13) {
    X_boston[[paste0(variable_names[i], "_x_", variable_names[j])]] <- X_boston[[i]] * X_boston[[j]]
  }
}

# Combine input and output
X_boston$Outcome <- y_boston

# Data Split
set.seed(42)  
train_indices <- sample(1:nrow(X_boston), size = 0.75 * nrow(X_boston))

# training and testing
train_input <- X_boston[train_indices, -ncol(X_boston)]  # Input variables for training
train_output <- X_boston[train_indices, ncol(X_boston)]  # Output variable for training
test_input <- X_boston[-train_indices, -ncol(X_boston)]  # Input variables for testing
test_output <- X_boston[-train_indices, ncol(X_boston)]  # Output variable for testing

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
  

  
  
  
  
  
  
  
  
  

