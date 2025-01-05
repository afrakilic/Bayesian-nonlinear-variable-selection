#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                                OZONE DATA                                   #####
#########################################################################################


#OZONE DATA 

library(mlbench)
data(Ozone)

Ozone <- Ozone[,-c(1,2,3, 9)] #removing the first three variables (month, day, week)
y_ozone <- Ozone$V4 #the outcome variable
Ozone <- Ozone[,-1] #removing the outcome variable from the matrix

ozone_results <- bayesian_selection(data_original = Ozone, knots = 6, y= y_ozone, iteration = 4000)

ozone_results$`selected model`
ozone_results$`posterior probability_selected`

#########################################################################################