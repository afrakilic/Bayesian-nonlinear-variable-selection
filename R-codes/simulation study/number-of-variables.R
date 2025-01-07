
######################################################################################################
#NUMBER OF VARIABLES (J)
######################################################################################################


n_var_values <- c(6, 10, 15, 20, 30)

results_n_var <- list()

# Run Bayesian selection for each n_var value
for (n_var in n_var_values) {
  # Create a unique name for each n_var value
  key <- paste0("J", n_var)
  
  # Store results in the list
  results_n_var[[key]] <- run_bayesian_selection( n_var = n_var)
}

######################################################################################################
#proportion of correct selection & posterior probability of the true model
J <- cbind(c(length(mj6[,1][mj6[,1]==1]), length(mj10[,1][mj10[,1]==1]),
             length(mj15[,1][mj15[,1]==1]),length(mj20[,1][mj20[,1]==1]),
             length(mj30[,1][mj30[,1]==1])), #proportion of correct selection
           
           c(mean(mj6[,2]), mean(mj10[,2]), mean(mj15[,2]), mean(mj20[,2]), mean(mj30[,2])),#posterior prob
           
           c(mj6_CPU, mj10_CPU, mj15_CPU, mj20_CPU, mj30_CPU)/mj10_CPU) #CPU  

colnames(J) = c("Correct Selection Proportion", "Posterior Pr. of the true model", "CPU")
rownames(J) = c("J=6", "J=10", "J=15", "J=20", "J=30")

J

######################################################################################################
multiplicity <- matrix(NA, ncol = 3, nrow = 5)
mj_results <- list(mj6, mj10, mj15, mj20, mj30)
for (i in 1:5){
  no_var <- matrix(NA, nrow = 100, ncol = 3)
  for (j in 1:100){
    result <- mj_results[[i]]
    no_var[j,] <- c(length(result[j,-c(1:3)][result[j,-c(1:3)]==1]), 
                    length(result[j,-c(1:3)][result[j,-c(1:3)]==2]),
                    length(result[j,-c(1:3)][result[j,-c(1:3)]==0]))
  }
  multiplicity[i,] <- colMeans(no_var)
}


#plot
plot(c(6, 10, 15, 20, 30), multiplicity[,1], type = "o", col = "darkblue",  ylim = c(0,27),#gamma=1
     main="Number of identified effect types across J", xlab="J", ylab="number of identified effect type") 
lines(c(6, 10, 15, 20, 30), multiplicity[,2], type = "o") #gamma=2
lines(c(6, 10, 15, 20, 30), multiplicity[,3], type = "s", lty = 2, pch = 15) #gamma=0
points(c(6, 10, 15, 20, 30), multiplicity[,3], pch = 18) #point for gamma=0
legend(5, 28, legend=c(expression(paste(gamma, "=1")), expression(paste(gamma, "=2")), expression(paste(gamma, "=0"))),
       col = c("darkblue", "black", "black"), lty=c(1,1,2), cex=1 )
