
#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                         BOSTON HOUSING DATA                                 #####
#########################################################################################


#data prep 
library(spdep)
data(boston)

#X Matrix
X_boston <- as.data.frame(cbind("CRIM" = boston.c$CRIM,
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


#Vector of responses 
y_boston <- (boston.c$CMEDV)

boston_results <- bayesian_selection(data_original = X_boston, y=y_boston, knots=6, iteration = 2000)

boston_results$`selected model`

bostonn <- boston_results
#########################################################################################


