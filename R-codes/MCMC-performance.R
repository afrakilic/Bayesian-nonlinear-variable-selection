#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######             MCMC MODEL SEARCH METHOD PERFORMANCE SIMULATION                  #####
#########################################################################################


#libraries 
library(mgcv)
library(data.table)
library(MCMCprecision)
library(lubridate)
library(progress)
  
#########################################################################################

#simulation function 
bayesian_selection<- function(beta = 0.8, n = 100, 
                              sigma2 = 0.1, n_var = 10, 
                              knots= 4,
                              iteration=1000,
                              gamma_prior = c(rep(0, n_var)),
                              prior_p=c(1/3, 1/3, 1/3)){
  
  #penalty is calculated as max edf - df(lm) for each variable swhere max edf = k-1 and df(lm)=1 
  penalty=knots-2
  #Data Generation
  #independent predictor variables 
  variables <- matrix(NA, nrow=n, ncol=n_var)
  colnames(variables) <-c(paste0("x", 1:n_var))
  for (i in 1:n_var){
    variables[,i] = rnorm(n)
  }
  #outcome variable
  error <- rnorm(n,sd=sqrt(sigma2)); y=error
  
  #randomly selecting the nonlinears and the zero relationships 
  nonlinear_true <- c("x1", "x2")
  linear_true <- c("x3", "x4")
  zero_true <- setdiff(colnames(variables), c(nonlinear_true, linear_true))
  
  # relationship definitions and creating the true model
  #nonlinear effects
  nonlinear_indices <- match(nonlinear_true, colnames(variables))
  for(i in 1:length(nonlinear_indices)){y=y+beta*exp(variables[, nonlinear_indices][,i])}
  
  #zero effects
  zero_indices <- match(zero_true, colnames(variables))
  for(i in 1:length(zero_indices)){y=y+0*variables[,zero_indices][,i]}
  
  #linear effects
  linear_indices<- match(linear_true, colnames(variables))
  for(i in 1: length(linear_indices)){y=y+beta*variables[,linear_indices][,i]}
  
  y = y - mean(y) #not considering the intercept 
  
  #true model gamma sequence 
  true_gamma = c(rep(NA, n_var))
  true_gamma[nonlinear_indices]=2; true_gamma[zero_indices]=0; true_gamma[linear_indices]=1
  data_original<- as.data.frame(variables) #transforming the matrix into a dataframe
  
  gamma_update_k <-gamma_prior #initial gamma specification
  gamma_draws <- matrix(NA, nrow= iteration, ncol = n_var) #matrix for gamma draws
  ps<- matrix(NA, nrow = iteration+1, ncol = 3) #matrix for p_draws
  
  pb<-progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = iteration,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar
  for(s in 1:iteration){
    pb$tick()
    for(k in 1:n_var){
      data = data_original
      gamma_update_k1 <- gamma_update_k[-c(k)]
      
      a <- data[, k] #the variable of interest 
      data <- data[, -k] #the remaining variables 
      
      # for linear effects
      linears <- c()
      if(length(gamma_update_k1[gamma_update_k1 == 1]) != 0){
        
        for(i in 1:ncol(data[which(gamma_update_k1 == 1)])){
          linears <- c(linears, paste(c(colnames(data[which(gamma_update_k1 == 1)][i])), collapse= ""))
        } 
      }
      
      #for nonlinear effects
      non_linears <- c()
      if(length(gamma_update_k1[gamma_update_k1 == 2]) != 0) {
        for(i in 1:ncol(data[which(gamma_update_k1 == 2)])){
          non_linears <- c(non_linears, paste(c('s(', colnames(data[which(gamma_update_k1 == 2)][i]), ',k=',knots,')'), collapse= ""))
        } 
      }
      
      vars <- c(linears, non_linears) #vector for non-zero effects to include to the models below
      
      #model fitting
      if(length(vars) != 0){ #remaining variables contain non-zero effect
        M1 <- gam(as.formula(paste('y', '~', paste(vars, collapse =  "+"))), data = data)
        M2 <- gam(as.formula(paste('y',  '~', 'a +', paste(vars, collapse =  "+"))), data = data) #1
        M3 <- gam(as.formula(paste('y',  '~', 's(a, k=',knots, ') +', paste(vars, collapse =  "+"))), data = data) #2
      } else{ #remaining variables do not contain non-zero effect
        M1 <- gam(y ~ 1, data = data)
        M2 <- gam(y ~ a, data = data) #1
        M3 <- gam(y ~ s(a, k=knots), data = data)
      }
      
      #BIC scores 
      bic_M1 <- (-2) * head(logLik(M1)) +  attr(logLik(M1), "df")* log(n)
      bic_M2 <- (-2) * head(logLik(M2)) +  attr(logLik(M2), "df")* log(n)
      bic_M3 <- (-2) * head(logLik(M3)) +  (attr(logLik(M2), "df") + penalty)* log(n) #penalty depends on the #of knots
      
      
      #Bayes Factors 
      BF11 <- exp((bic_M1 - bic_M1) /2) #null against the null
      BF21_ <- exp((bic_M1 - bic_M2) /2) #linear against the null
      BF31_ <- exp((bic_M1 - bic_M3) /2) #nonlinear against the null 
      
      
      #infinity BFs
      if(BF21_ == "-Inf"){BF21 = -1e5} else if(BF21_== "Inf") {BF21 = 1e5} else{BF21=BF21_}
      if(BF31_ == "-Inf"){BF31 = -1e5} else if(BF31_ == "Inf") {BF31 = 1e5} else{BF31=BF31_}
      
      #Posterior Probabilities
      zero <- (prior_p[1] * BF11) / sum(prior_p[1] * BF11, prior_p[2]*BF21, prior_p[3]* BF31)
      one <- (prior_p[2] * BF21) / sum(prior_p[1] * BF11, prior_p[2]*BF21, prior_p[3]* BF31)
      two <- (prior_p[3] * BF31) / sum(prior_p[1] * BF11, prior_p[2]*BF21, prior_p[3]* BF31)
      
      
      
      #sampling the effect type of the variable of interest 
      gamma_update_k[k] <- sample(c(0,1,2), size=1, prob = c(zero, one, two)) 
    }
    
    #multiplicity correction for the next chain 
    prior_p <- rdirichlet(1, c(1+length(gamma_update_k[which(gamma_update_k == 0)]), #alpha1 = 1
                               1+length(gamma_update_k[which(gamma_update_k == 1)]),  #alpha2 = 1
                               1+length(gamma_update_k[which(gamma_update_k == 2)])))
    ps[s+1, ] <- prior_p
    gamma_draws[s,] <- gamma_update_k
  }
  #posterior probability calculation of the true model 
  posterior <- data.table(gamma_draws)
  frequency <- as.matrix(posterior[,list(posterior=.N),by=names(posterior)][order(posterior,decreasing=T)])
  t=apply(frequency[,-(n_var+1)], 1, function(x) return(all(x == true_gamma))) #checking any true model among the draws
  pp_t<- if (any(t)==FALSE) 0 else as.numeric(frequency[which(t),(n_var+1)])/iteration ##posterior prob calculation of the true
  pp_s <- as.numeric(frequency[1,(n_var+1)])/iteration ##posterior prob calculation of the selected
  results <- list("is true" = if(identical(as.vector(frequency[1,1:n_var]),true_gamma) == TRUE) 1 else 0,
                  "posterior probability_true" = pp_t,
                  "posterior probability_selected" = pp_s,
                  "selected model" = as.vector(frequency[1,1:n_var]),
                  "true model" = true_gamma,
                  "gamma draws" = gamma_draws,
                  "p draws" = ps)
  return(results)
}



#########################################################################################

# CONVERGENCE OF GIBSS SAMPLER CHAIN 

#J=6
chain_6 <- bayesian_selection_m(iteration = 10000, n_var = 6) #fit

chain_6_pp <- c(rep(NA, 100))#matrix to store the posterior probabilities at every 100 iteration
for (i in 1:100){
  gamma_draws <- chain_6$`gamma draws`
  posterior <- data.table(gamma_draws[1:(i*100),])
  frequency <- as.matrix(posterior[,list(posterior=.N),by=names(posterior)][order(posterior,decreasing=T)])
  pp_s <- as.numeric(frequency[1,(6+1)])/(i*100) ##posterior prob calculation of the selected
  chain_6_pp[i] = pp_s
}

#plot
plot(c(100*c(1:100)), chain_6_pp, type="l", ylim=c(0, 1), 
     ylab = "Posterior probability of the selected model", xlab = "T", main = "J=6")

######################################################################################################
#J=20
chain_20 <- bayesian_selection_m(iteration = 10000, n_var = 20) #fit

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


# Run the loop in parallel
detectCores()
#32
registerDoParallel(30)
trials=100


######################################################################################################
                                    #SAMPLE SIZES (n)
######################################################################################################



xresults_n50 <- foreach(icount(trials),
                        .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                        .combine = rbind) %dopar%{
                          # set the specific values of the parameters 
                          a = bayesian_selection(beta=.5, n=50)
                          b = bayesian_selection(beta=2, n=50)
                          
                          c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                            a$`selected model`,
                            b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                            b$`selected model`)
                        }

######################################################################################################
xresults_n100 <- foreach(icount(trials),
                         .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                         .combine = rbind) %dopar%{
                           # set the specific values of the parameters 
                           a = bayesian_selection(beta=.5, n=100)
                           b = bayesian_selection(beta=2, n=100)
                           
                           c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                             a$`selected model`,
                             b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                             b$`selected model`)
                         }
######################################################################################################
xresults_n250 <- foreach(icount(trials),
                         .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                         .combine = rbind) %dopar%{
                           # set the specific values of the parameters 
                           a = bayesian_selection(beta=.5, n=250)
                           b = bayesian_selection(beta=2, n=250)
                           
                           c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                             a$`selected model`,
                             b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                             b$`selected model`)
                         }
######################################################################################################
xresults_n500 <- foreach(icount(trials),
                         .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                         .combine = rbind) %dopar%{
                           # set the specific values of the parameters 
                           a = bayesian_selection(beta=.5, n=500)
                           b = bayesian_selection(beta=2, n=500)
                           
                           c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                             a$`selected model`,
                             b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                             b$`selected model`)
                         }


######################################################################################################
#proportion of correct selection  

#beta=.5
cs_n_.5 <- c(length(n50[,1][n50[,1]==1]), 
             length(n100[,1][n100[,1]==1]), 
             length(n250[,1][n250[,1]==1]),
             length(n500[,1][n500[,1]==1]))

#beta=2
cs_n_2 <- c(length(n50[,14][n50[,14]==1]), 
            length(n100[,14][n100[,14]==1]), 
            length(n250[,14][n250[,14]==1]),
            length(n500[,14][n500[,14]==1]))

#plot
plot(c( 50, 100, 250, 500), cs_n_.5/100, type="l", ylim=c(0, 1), 
     ylab = "Proportion of correct selection", xlab = "n")
legend(400, 0.5, legend=c(expression(paste(beta, "=.5")), expression(paste(beta, "=2"))),
       lty=1:2, cex=1)
lines(c( 50, 100, 250, 500), cs_n_2/100, type="l", lty=2)


######################################################################################################
#posterior probability of the true model 

#beta=.5
pp_n_.5 <-c(mean(n50[,2]),
            mean(n100[,2]),
            mean(n250[,2]),
            mean(n500[,2])) #
#beta=2
pp_n_2 <- c(mean(n50[,15]),
            mean(n100[,15]),
            mean(n250[,15]),
            mean(n500[,15]))

#plot
plot(c( 50, 100, 250, 500), pp_n_.5, type="l", ylim=c(0, 1), 
     ylab = "Posterior Probability of the true model", xlab = expression(beta))
legend(400, 0.5, legend=c(expression(paste(beta, "=.5")), expression(paste(beta, "=2"))),
       lty=1:2, cex=1)
lines(c( 50, 100, 250, 500), pp_n_2, type="l", lty=2)



######################################################################################################
                                    #EFFECT SIZES (BETA)
######################################################################################################



#beta=0.1
xresults_beta_.1 <- foreach(icount(trials),
                            .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                            .combine = rbind) %dopar%{
                              # set the specific values of the parameters 
                              a = bayesian_selection(beta=.1, n=100)
                              b = bayesian_selection(beta=.1, n=250)
                              
                              c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                                a$`selected model`,
                                b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                                b$`selected model`)
                            }


######################################################################################################
#beta=0.5
xresults_beta_.5 <- foreach(icount(trials),
                            .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                            .combine = rbind) %dopar%{
                              # set the specific values of the parameters 
                              a = bayesian_selection(beta=.5, n=100)
                              b = bayesian_selection(beta=.5, n=250)
                              
                              c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                                a$`selected model`,
                                b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                                b$`selected model`)
                            }


######################################################################################################
#beta=0.1
xresults_beta_1 <- foreach(icount(trials),
                           .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                           .combine = rbind) %dopar%{
                             # set the specific values of the parameters 
                             a = bayesian_selection(beta=1, n=100)
                             b = bayesian_selection(beta=1, n=250)
                             
                             c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                               a$`selected model`,
                               b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                               b$`selected model`)
                           }


######################################################################################################
#beta=0.1
xresults_beta_1.5 <- foreach(icount(trials),
                             .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                             .combine = rbind) %dopar%{
                               # set the specific values of the parameters 
                               a = bayesian_selection(beta=1.5, n=100)
                               b = bayesian_selection(beta=1.5, n=250)
                               
                               c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                                 a$`selected model`,
                                 b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                                 b$`selected model`)
                             }



######################################################################################################
#beta=0.1
xresults_beta_2 <- foreach(icount(trials),
                           .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                           .combine = rbind) %dopar%{
                             # set the specific values of the parameters 
                             a = bayesian_selection(beta=2, n=100)
                             b = bayesian_selection(beta=2, n=250)
                             
                             c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                               a$`selected model`,
                               b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                               b$`selected model`)
                           }


######################################################################################################
xresults_beta_2.5 <- foreach(icount(trials),
                             .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                             .combine = rbind) %dopar%{
                               # set the specific values of the parameters 
                               a = bayesian_selection(beta=2.5, n=100)
                               b = bayesian_selection(beta=2.5, n=250)
                               
                               c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                                 a$`selected model`,
                                 b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                                 b$`selected model`)
                             }

######################################################################################################
xresults_beta_3 <- foreach(icount(trials),
                           .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                           .combine = rbind) %dopar%{
                             # set the specific values of the parameters 
                             a = bayesian_selection(beta=3, n=100)
                             b = bayesian_selection(beta=3, n=250)
                             
                             c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                               a$`selected model`,
                               b$`is true`, b$`posterior probability_true`,b$`posterior probability_selected`,
                               b$`selected model`)
                           }


######################################################################################################
#proportion of correct selection

#n=100
cs_beta_100 <- c(length(beta.1[,1][beta.1[,1]==1]), 
                 length(beta.5[,1][beta.5[,1]==1]), 
                 length(beta1[,1][beta1[,1]==1]),
                 length(beta1.5[,1][beta1.5[,1]==1]),
                 length(beta2[,1][beta2[,1]==1]),
                 length(beta2.5[,1][beta2.5[,1]==1]), 
                 length(beta3[,1][beta3[,1]==1]))

#n=250
cs_beta_250 <- c(length(beta.1[,14][beta.1[,14]==1]),
                 length(beta.5[,14][beta.5[,14]==1]),
                 length(beta1[,14][beta1[,14]==1]),
                 length(beta1.5[,14][beta1.5[,14]==1]),
                 length(beta2[,14][beta2[,14]==1]),
                 length(beta2.5[,14][beta2.5[,14]==1]),
                 length(beta3[,14][beta3[,14]==1]))



#plot
plot(c( .1,.5,1, 1.5, 2, 2.5, 3), cs_beta_100/100, type="l", ylim=c(0, 1), 
     ylab = "Proportion of correct selection", xlab = expression(beta))
legend(2, 0.3, legend=c("n=100", "n=250"),
       lty=1:2, cex=1)
lines(c( .1,.5,1, 1.5, 2, 2.5, 3), cs_beta_250/100, type="l", lty=2)


######################################################################################################
#posterior probability of the true model 

#n=100
pp_beta_100 <- c(mean(beta.1[,2]),
                 mean(beta.5[,2]), 
                 mean(beta1[,2]), 
                 mean(beta1.5[,2]),
                 mean(beta2[,2]),
                 mean(beta2.5[,2]), 
                 mean(beta3[,2]))

#n=250
pp_beta_250 <- c(mean(beta.1[,15]),
                 mean(beta.5[,15]),
                 mean(beta1[,15]),
                 mean(beta1.5[,15]),
                 mean(beta2[,15]),
                 mean(beta2.5[,15]),
                 mean(beta3[,15]))

#plot
plot(c( .1,.5,1, 1.5, 2, 2.5, 3), pp_beta_100, type="l", ylim=c(0, 1), 
     ylab = "Posterior Probability of the true model", xlab = expression(beta))
legend(2, 0.3, legend=c("n=100", "n=250"),
       lty=1:2, cex=1)
lines(c( .1,.5,1, 1.5, 2, 2.5, 3), pp_beta_250, type="l", lty=2)


######################################################################################################
                                  #NUMBER OF VARIABLES (J)
######################################################################################################




start <- proc.time()
mj6 <- foreach(icount(100),
               .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
               .combine = rbind) %dopar%{
                 # set the specific values of the parameters 
                 a = bayesian_selection(n_var=6)
                 
                 
                 c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                   a$`selected model`)
               }

mj6_CPU <- proc.time()-start

######################################################################################################

start <- proc.time()
mj10 <- foreach(icount(100),
                .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                .combine = rbind) %dopar%{
                  # set the specific values of the parameters 
                  a = bayesian_selection(n_var=10)
                  
                  
                  c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                    a$`selected model`)
                }

mj10_CPU <- proc.time()-start


######################################################################################################
start <- proc.time()
mj20 <- foreach(icount(100),
                .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                .combine = rbind) %dopar%{
                  # set the specific values of the parameters 
                  a = bayesian_selection(n_var=20)
                  
                  
                  c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                    a$`selected model`)
                }

mj20_CPU <- proc.time()-start

######################################################################################################

start <- proc.time()
mj15 <- foreach(icount(100),
                .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                .combine = rbind) %dopar%{
                  # set the specific values of the parameters 
                  a = bayesian_selection(n_var=15)
                  
                  
                  c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                    a$`selected model`)
                }

mj15_CPU <- proc.time()-start

######################################################################################################

start <- proc.time()
mj30 <- foreach(icount(100),
                .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                .combine = rbind) %dopar%{
                  # set the specific values of the parameters 
                  a = bayesian_selection(n_var=30)
                  
                  
                  c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                    a$`selected model`)
                }

mj30_CPU <- proc.time()-start

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


######################################################################################################
                                      #NUMBER OF BASIS FUNCTIONS (k)
######################################################################################################

start <- proc.time()
knots4 <- foreach(icount(100),
                  .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                  .combine = rbind) %dopar%{
                    # set the specific values of the parameters 
                    a = bayesian_selection(nots=4)
                    
                    
                    c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                      a$`selected model`)
                  }

write.csv(nots4, "nots4.csv")
knots4_CPU <- proc.time()-start

######################################################################################################

start <- proc.time()
knots6 <- foreach(icount(100),
                  .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                  .combine = rbind) %dopar%{
                    # set the specific values of the parameters 
                    a = bayesian_selection(nots=6)
                    
                    c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                      a$`selected model`)
                  }

write.csv(nots6, "nots6.csv")
knots6_CPU <- proc.time()-start

######################################################################################################

start <- proc.time()
knots8 <- foreach(icount(100),
                  .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                  .combine = rbind) %dopar%{
                    # set the specific values of the parameters 
                    a = bayesian_selection(nots=8)
                    
                    
                    c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                      a$`selected model`)
                  }

write.csv(nots8, "nots8.csv")
knots8_CPU <- proc.time()-start

######################################################################################################

start <- proc.time()
knots10 <- foreach(icount(100),
                   .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                   .combine = rbind) %dopar%{
                     # set the specific values of the parameters 
                     a = bayesian_selection(nots=10)
                     
                     
                     c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                       a$`selected model`)
                   }

knots10_CPU <- proc.time()-start

######################################################################################################

start <- proc.time()
knots12 <- foreach(icount(100),
                   .packages = c("mgcv", "data.table", "MCMCprecision", "lubridate", "progress"),
                   .combine = rbind) %dopar%{
                     # set the specific values of the parameters 
                     a = bayesian_selection(nots=12)
                     
                     
                     c(a$`is true`, a$`posterior probability_true`,a$`posterior probability_selected`,
                       a$`selected model`)
                   }

knots12_CPU <- proc.time()-start


######################################################################################################

#proportion of correct selection & posterior probability of the true model
knots <- cbind(c(length(knots4[,1][knots4[,1]==1]), length(knots6[,1][knots6[,1]==1]),
                 length(knots8[,1][knots8[,1]==1]),length(knots10[,1][knots10[,1]==1]),
                 length(knots12[,1][knots12[,1]==1])),
               
               c(mean(knots4[,2]), mean(knots6[,2]),
                 mean(knots8[,2]), mean(knots10[,2]), mean(knots12[,2])),
               c(mj10_CPU, knots6_CPU, knots8_CPU, knots10_CPU, knots12_CPU)/mj10_CPU) #CPU  

colnames(knots) = c("Correct Selection Proportion", "Posterior Pr. of the true model", "CPU")

rownames(knots) = c("k=4", "k=6", "k=8", "k=10", "k=12")

knots



######################################################################################################

######################################################################################################





