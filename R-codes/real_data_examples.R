#REAL DATA EXAMPLES 

#libraries 
library(ISLR)
library(mgcv)


#COLLEGE DATA 
data_original <- as.data.frame(College[,-c(,3,4)]) #removing  Enroll, Accept due to potential causal relationship with number of applications.
y <- data_original$Apps

data_original <- data_original[,-1]
  n=777
  n_var = 14
  nots= 4
  penalty=nots-2
  iteration=2000
  gamma_prior = c(rep(0, n_var))
  prior_p=c(1/3, 1/3, 1/3)

  
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
          non_linears <- c(non_linears, paste(c('s(', colnames(data[which(gamma_update_k1 == 2)][i]), ',k=',nots,')'), collapse= ""))
        } 
      }
      
      vars <- c(linears, non_linears) #vector for non-zero effects to include to the models below
      
      #model fitting
      if(length(vars) != 0){ #remaining variables contain non-zero effect
        M1 <- gam(as.formula(paste('y', '~', paste(vars, collapse =  "+"))), data = data)
        M2 <- gam(as.formula(paste('y',  '~', 'a +', paste(vars, collapse =  "+"))), data = data) #1
        M3 <- gam(as.formula(paste('y',  '~', 's(a, k=',nots, ') +', paste(vars, collapse =  "+"))), data = data) #2
      } else{ #remaining variables do not contain non-zero effect
        M1 <- gam(y ~ 1, data = data)
        M2 <- gam(y ~ a, data = data) #1
        M3 <- gam(y ~ s(a, k=nots), data = data)
      }
      
      #BIC scores 
      bic_M1 <- (-2) * head(logLik(M1)) +  attr(logLik(M1), "df")* log(n)
      bic_M2 <- (-2) * head(logLik(M2)) +  attr(logLik(M2), "df")* log(n)
      bic_M3 <- (-2) * head(logLik(M3)) +  (attr(logLik(M2), "df") + penalty)* log(n) #penalty depends on the #of nots
      
      
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

  plot(data_original$Top10perc, y)
  
  plot(College$Private,y)
College$Apps  
