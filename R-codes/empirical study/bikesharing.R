
#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                         BIKE-SHARING DATA (DAILY)                            #####
#########################################################################################


source("R-codes/bayesian_selection.R")  # Ensure the script is in the working directory
source("R-codes/utils.R")
library(glmnet)


prelim.data=read.table("/Users/hakilic/Downloads/day.csv",header = TRUE,sep = ",")

# remove date and record index since these variables are just indexing variables
# remove casual and registered since they define the outcome variable of interest
prelim.data=prelim.data[ ,-c(1,2,14,15)]

# Changing normalised versions to actual numbers based on the read me file
prelim.data$temp <- prelim.data$temp*41
prelim.data$atemp <- prelim.data$atemp*50
prelim.data$windspeed <- prelim.data$windspeed*67
prelim.data$hum <- prelim.data$hum*100

prelim.data$yr <- factor(prelim.data$yr, levels = c(0, 1), labels = c("2011", "2012"))
prelim.data$weathersit <- factor(format(prelim.data$weathersit, format="%A"),
                                 levels = c("1", "2","3") ,
                                 labels = c("Good","Moderate","Bad"))
prelim.data$holiday <- factor(format(prelim.data$holiday, format="%A"),
                              levels = c("0", "1") , labels = c("NotHoliDay","Holiday"))
prelim.data$season <- factor(format(prelim.data$season, format="%A"),
                             levels = c("1", "2","3","4") , labels = c("Spring","Summer","Fall","Winter"))
prelim.data$mnth <- factor(prelim.data$mnth)
prelim.data$mnth <- relevel(prelim.data$mnth,ref=6)
prelim.data$weekday <- factor(prelim.data$weekday)
prelim.data$weekday <- relevel(prelim.data$weekday,ref=3)
prelim.data=subset(prelim.data,select = -c(workingday))
prelim.data$cnt=(prelim.data$cnt)^(1/2)

dummy_variables <- as.data.frame(model.matrix(~ yr + weathersit + holiday + season + mnth + weekday - 1, data = prelim.data))

# centering the Y. Assuming Y is always the last column of prelim.data
y.val.0 <- mean(prelim.data[ , ncol(prelim.data)])
prelim.data[ , ncol(prelim.data)]=prelim.data[ , ncol(prelim.data)]-mean(prelim.data[ , ncol(prelim.data)])

X <- cbind(dummy_variables, prelim.data[,7:11])

# Data Split
set.seed(42)  
train_indices <- sample(1:nrow(X), size = 0.75 * nrow(X))

# training and testing
train_input <- X[train_indices, -ncol(X)]  # Input variables for training
train_output <- X[train_indices, ncol(X)]  # Output variable for training
test_input <- X[-train_indices, -ncol(X)]  # Input variables for testing
test_output <- X[-train_indices, ncol(X)]  # Output variable for testing

# Training
bikesharing_results <- bayesian_selection(X = train_input, y=train_output, knots=10, iteration = 2000)

selected_model <- bikesharing_results$`selected model`

results <- evaluate_model(
  selected_model =selected_model,
  train_input = train_input,
  train_output = train_output,
  test_input = test_input,
  test_output = test_output,
  knots = 10
)
