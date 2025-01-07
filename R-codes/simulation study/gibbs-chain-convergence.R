#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######             MCMC MODEL SEARCH METHOD PERFORMANCE SIMULATION                  #####
#########################################################################################

source("R-codes/bayesian_selection.R")
source("R-codes/utils.R")
#libraries 
library(mgcv)
library(data.table)
library(MCMCprecision)
library(lubridate)
library(progress)
  
#########################################################################################


#########################################################################################

# CONVERGENCE OF GIBSS SAMPLER CHAIN 

#J=6

data <- generate_true_model(n_var = 6)
chain_6 <- bayesian_selection(X = data$data, y= data$response, iteration = 10000) #fit


chain_6_pp <- c(rep(NA, 100))#matrix to store the posterior probabilities at every 100 iteration
for (i in 1:100){
  gamma_draws <- chain_6$`gamma draws`
  posterior <- data.table(gamma_draws[1:(i*100),])
  frequency <- as.matrix(posterior[,list(posterior=.N),by=names(posterior)][order(posterior,decreasing=T)])
  pp_s <- as.numeric(frequency[1,(6+1)])/(i*100) ##posterior prob calculation of the selected
  chain_6_pp[i] =  as.numeric(frequency[1,(6+1)])/(i*100)
}

#plot
plot(c(100*c(1:100)), chain_6_pp, type="l", ylim=c(0, 1), 
     ylab = "Posterior probability of the selected model", xlab = "T", main = "J=6")


######################################################################################################
#J=20
data <- generate_true_model(n_var = 20)
chain_20 <- bayesian_selection(X = data$data, y= data$response, iteration = 10000) #fit

chain_20_pp <- c(rep(NA, 100))#matrix to store the posterior probabilities at every 100 iteration
for (i in 1:100){
  gamma_draws <- chain_20$`gamma draws`
  posterior <- data.table(gamma_draws[1:(i*100),])
  frequency <- as.matrix(posterior[,list(posterior=.N),by=names(posterior)][order(posterior,decreasing=T)])
  pp_s <- as.numeric(frequency[1,(20+1)])/(i*100) ##posterior prob calculation of the selected
  chain_20_pp[i] = pp_s
}

#plot
plot(c(100*c(1:100)), chain_20_pp, type="l", ylim=c(0, 1), 
     ylab = "Posterior probability of the selected model", xlab = "T", main = "J=20")


#########################################################################################







