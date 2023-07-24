#### Title: Bayesian Variable Selection for linear and Nonlinear Model
#### Author: Afra Kilic
#### Created: August, 2023

#########################################################################################
######                 BAYES FACTOR CONSISTENCY SIMULATION                          #####
#########################################################################################


#libraries 
library(mgcv)

##########################################################################################
                              #NULL MODEL TEST
##########################################################################################


null_model_test <- function(n = 100, sigma2 = 0.1){ #default n=100 and sigma2=0.1
  
  output_mat_null <- matrix(NA, nrow = 500, ncol = 2)
  colnames(output_mat_null) <- c("B01", "BF02")
  
  for (j in 1:500){ #for 500 different datasets 
    data <- as.data.frame(cbind("y" = rnorm(n, 1, sigma2),#data generation 
                                "x1" = rnorm(n, 0, 0.5)))
    
    #model fit
    null_model <- lm(y ~ 1, data = data)
    linear_model <- lm(y ~  1 + .,data = data)
    non_linear_model <- gam(y~ 1 + s(x1, k=4), data = data)
    
    #BIC calculation 
    bic_null <- (-2) * head(logLik(null_model)) +  attr(logLik(null_model), "df")* log(n) #H0 
    bic_lin <- (-2) * head(logLik(linear_model)) +  attr(logLik(linear_model), "df")* log(n) #H1
    bic_nonl <- (-2) * head(logLik(non_linear_model)) +  (attr(logLik(linear_model), "df") + 2) * log(n) #H2
    
    
    B01 <- log(exp((bic_lin - bic_null) /2)) #null against the linear
    BF02 <- log(exp((bic_nonl - bic_null) /2)) #null against the nonlinear model 
    
    output_mat_null[j, 1] <- B01
    output_mat_null[j, 2] <- BF02
  }
  
  return(output_mat_null)
}  


##########################################################################################

#Sample Sizes

n_sizes <- c(50, 100, 500, 1000) #sample sizes 
null_samp <- matrix(NA, nrow = 500, ncol = 8)  #matrix for the results 
colnames(null_samp) <- c("B01_50","BF02_50", "B01_100", "BF02_100", 
                         "B01_500", "BF02_500", "B01_1000", "BF02_1000")
for (i in 1:4){
  null_samp[,c((i + i-1), (2*i))] <- null_model_test(n=n_sizes[i])
}


#PLots 
plot(density(null_samp[,1]), 
     xlim = c(min(null_samp), max(null_samp)), ylim=c(0, 1.7),
     main = "log(BF01) for different sample sizes", col = "red",
     xlab = "log(BF01)")
lines(density(null_samp[,3]), col = "darkblue")
lines(density(null_samp[,5]), col = "magenta")
lines(density(null_samp[,7]))
legend("topright", col = c("red", "darkblue", "magenta", "black"), cex=.80,
       lty=c(1,2,1,1),lwd=c(2,2,1,1), inset=.01,
       legend = c("n = 50", "n = 100", "n = 500", "n = 1000"))


plot(density(null_samp[,2]),
     xlim = c(min(null_samp), max(null_samp)), ylim=c(0, 1.7),
     main = "log(BF02) for different sample sizes", col = "red",
     xlab = "log(BF02)")
lines(density(null_samp[,4]), col = "darkblue")
lines(density(null_samp[,6]), col = "magenta")
lines(density(null_samp[,8]))
legend("topright", col = c("red", "darkblue", "magenta", "black"), cex=.80,
       lty=c(1,2,1,1),lwd=c(2,2,1,1), inset=.01,
       legend = c("n = 50", "n = 100", "n = 500", "n = 1000"))



##########################################################################################

#Sigma 

sigma2 <- c(seq(from = 0.01,  to = 1, length.out = 10)) #sigma2 sequence

null_sigma <- matrix(NA, nrow = 500, ncol = 20) #matrix for the results 
colnames(null_sigma) <- c("B01, s2=0.01", "BF02, s2=0.01",
                          "B01, s2=0.12", "BF02, s2=0.12",
                          "B01, s2=0.23", "BF02, s2=0.23",
                          "B01, s2=0.34", "BF02, s2=0.34",
                          "B01, s2=0.45", "BF02, s2=0.45",
                          "B01, s2=0.56", "BF02, s2=0.56",
                          "B01, s2=0.67", "BF02, s2=0.67",
                          "B01, s2=0.78", "BF02, s2=0.78",
                          "B01, s2=0.89", "BF02, s2=0.89",
                          "B01, s2=1.00", "BF02, s2=1.00")
for (i in 1:10){
  null_sigma[, c((i + i-1), (2*i))] <- null_model_test(sigma2 = sigma2[i])
}

plot(x = sigma2, y=colMeans(null_sigma[,c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19)]),
     type = "l", col="red", xlab = "sigma2", ylab = "log(B01)",  
     main = "log(BF01) for different sigma2")

plot(x = sigma2, y=colMeans(null_sigma[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)]),
     type = "l", col="blue", xlab = "sigma2", ylab = "log(BF02)",
     main = "log(BF02) for different sigma2")

##########################################################################################
                                  #LINEAR MODEL TEST
##########################################################################################

linear_model_test <- function(beta = 1.5, n = 100, sigma2 = 0.1){
  BF01s <- c(rep(NA, 500))
  for (j in 1:500){
    error <- rnorm(n,sd=sqrt(sigma2))
    x1 = rnorm(n, 0, 0.5)
    y <- beta*x1 + error
    data <- as.data.frame(cbind(y, x1))
    
    #model fit
    linear_model <- lm(y ~ 1 + x1, data = data) #H0 
    non_linear_model <- gam(y~ 1 + s(x1, k=4), data = data) #H1
    
    #BIC calculation 
    bic_lin <- (-2) * head(logLik(linear_model)) + 3* log(n) #H0
    bic_nonl <- (-2) * head(logLik(non_linear_model)) +  5 * log(n) #H1
    
    logLik(linear_model);logLik(non_linear_model)
    
    BF01 <- log(exp((bic_nonl - bic_lin) /2)) #BF_{linear against nonlinear} higher values indicates more evidence for the linear model
    BF01s[j] <- BF01}
  
  return(BF01s)
}

##########################################################################################

#Sample Size
linear_samp <- matrix(NA, nrow = 500, ncol = 4) #matrix for the results  
colnames(linear_samp) <- c("BF01_50","BF01_100", "BF01_500", "BF01_1000")

for (i in 1:4){
  linear_samp[,i] <- linear_model_test(n = n_sizes[i])
}

##########################################################################################

#Effect Size (beta)
betas <- seq(-3, 3, length.out = 10) #beta sequence

linear_beta <- matrix(NA, nrow=500, ncol=10) #matrix for the results 
for (i in 1:10){
  linear_beta[,i] <- linear_model_test(beta = betas[i])
}

plot(y=colMeans(linear_beta), x=betas, type="l", ylab="log(BF01)", xlab="beta", main="log(BF01) across Beta values")


##########################################################################################

#Error Rate (sigma)
sigma2 <- c(seq(from = 0.01,  to = 1, length.out = 10)) #sigma2

linear_sigma <- matrix(NA, nrow=500, ncol=10)
colnames(linear_sigma) <- as.character(sigma2)
for (i in 1:10){
  linear_sigma[,i] <- linear_model_test(sigma2 = sigma2[i])
}

plot(y=colMeans(linear_sigma), x=sigma2, type = "l", xlab="sigma2", ylab="log(BF01)", main="log(BF01) across sigma2 values")


##########################################################################################
                                  #NONLINEAR MODEL TEST
##########################################################################################
non_linear_test <- function(beta = 1.5, sigma2 = 0.1, n= 100){
  #default beta= 1, n=100 and sigma2=0.5
  output_mat_nonl <- matrix(NA, nrow = 500, ncol = 6)
  colnames(output_mat_nonl) <- c("BF01_lin", "BF01_log", "BF01_exp", "BF01_quad", "BF01_sine", "BF01_norm")
  for (j in 1:500){
    
    x1 = rnorm(n, 0, 0.5)
    error <- rnorm(n,sd=sqrt(sigma2)) #error terms
    y_lin <- beta*x1 + error  #linear 
    y_log <- beta*log(x1) + error #logarithmic
    y_exp <- beta*exp(x1) + error #exponential
    y_quad <- beta*x1^4 + error #quadratic
    y_sine <- beta*sin(x1/3) + error #sine
    y_norm <- beta*dnorm(x1)*x1 + error #a nonlinear model
    data <- as.data.frame(cbind(y_lin, y_log, y_exp, y_quad, y_sine, y_norm, x1))
    
    #################################################################################
    
    #Model fits
    
    #linear
    linear_model_lin <- lm(y_lin ~ 1 + x1, data=data)
    non_linear_model_lin <- gam(y_lin ~ 1 + s(x1, k= 4), data=data)
    non_linear_model_lin$edf
    #logarithmic
    linear_model_log <- lm(y_log ~ 1 + x1, data = data)
    non_linear_model_log <- gam(y_log ~ 1 + s(x1, k=4), data=data)
    
    #exponential
    linear_model_exp <- lm(y_exp ~ 1 + x1, data = data)
    non_linear_model_exp <- gam(y_exp ~ 1 + s(x1, k=4), data=data)
    
    #quadratic
    linear_model_quad <- lm(y_quad ~ 1 + x1, data = data)
    non_linear_model_quad <- gam(y_quad ~ 1 +  s(x1, k=4), data=data)
    
    #sine
    linear_model_sine <- lm(y_log ~ 1 + x1, data = data)
    non_linear_model_sine <- gam(y_log~ 1 + s(x1, k=4), data=data)
    
    #a-nonlinear
    linear_model_norm <- lm(y_log ~ 1 + x1, data = data)
    non_linear_model_norm <- gam(y_log~ 1 + s(x1, k=4), data=data)
    
    #################################################################################
    
    #BIC calculations
    
    #linear
    bic_lin <- (-2) * head(logLik(linear_model_lin)) + 3*log(n) #H0
    bic_nonl <- (-2) * head(logLik(non_linear_model_lin)) + 5*log(n) #H1
    BF01_lin <- log(exp((bic_nonl - bic_lin) /2))
    output_mat_nonl[j,1] <- BF01_lin
    
    #logarithmic
    bic_lin_log <- (-2) * head(logLik(linear_model_log)) +  3*log(n) #H0
    bic_nonl_log <- (-2) * head(logLik(non_linear_model_log)) + 5*log(n) #H1
    BF01_log <- log(exp((bic_nonl_log - bic_lin_log) /2))
    output_mat_nonl[j,2] <- BF01_log
    
    #exponential 
    bic_lin_exp <- (-2) * head(logLik(linear_model_exp)) +  3*log(n) #H0
    bic_nonl_exp <- (-2) * head(logLik(non_linear_model_exp)) +  5*log(n) #H1
    BF01_exp <- log(exp((bic_nonl_exp - bic_lin_exp) /2))
    output_mat_nonl[j,3] <- BF01_exp
    
    #quadratic
    bic_lin_quad <- (-2) * head(logLik(linear_model_quad)) +  3*log(n) #H0
    bic_nonl_quad <- (-2) * head(logLik(non_linear_model_quad)) + 5*log(n) #H1
    BF01_quad <- log(exp((bic_nonl_quad - bic_lin_quad) /2))
    output_mat_nonl[j,4] <- BF01_quad
    
    #sine
    bic_lin_sine <- (-2) * head(logLik(linear_model_sine)) + 3*log(n) #H0
    bic_nonl_sine <- (-2) * head(logLik(non_linear_model_sine)) + 5*log(n) #H1
    BF01_sine <- log(exp((bic_nonl_sine - bic_lin_sine) /2))
    output_mat_nonl[j,5] <- BF01_sine
    
    #a-nonlinear
    bic_lin_norm <- (-2) * head(logLik(linear_model_norm)) +3*log(n) #H0
    bic_nonl_norm <- (-2) * head(logLik(non_linear_model_norm)) + 5*log(n) #H1
    BF01_norm <- log(exp((bic_nonl_norm - bic_lin_norm) /2))
    output_mat_nonl[j,6] <- BF01_sine
    
  }
  return(output_mat_nonl)
}

##########################################################################################

#Sample Sizes 

B0_samp <- c() #vector for the results 
for (i in 1: 4){
  samp <- non_linear_test(n = n_sizes[i])
  B0_samp <- cbind(B0_samp, samp)
}

n_means<- colMeans(B0_samp)
plot(c(seq(from =50, to = 1000, length = 4)), n_means[c(1, 7, 13, 19)], ylim = c(min(n_means),max(n_means)),
     type = "l", col = "red", main = "log(BF01) for different sample sizes", 
     ylab = "log(BF01)", xlab = "Sample Size (n)") #linear
lines(c(seq(from =50, to = 1000, length = 4)), n_means[c(2, 8, 14, 20)], col = "darkblue") #logarithmic
lines(seq(from =50, to = 1000, length = 4), n_means[c(3, 9, 15, 21)], col = "magenta") #exponential
lines(seq(from =50, to = 1000, length = 4), n_means[c(4, 10, 16, 22)], col= "orange") #quadratic
lines(seq(from =50, to = 1000, length = 4), n_means[c(5, 11, 17, 23)], col = "pink") #sine
lines(seq(from =50, to = 1000, length = 4), n_means[c(6, 12, 18, 24)]) #a-nonlinear 
lines(seq(from =50, to = 1000, length = 4), y = c(rep(2.3, 4)),  col = "green", lwd = 2, lty = 2) #BF01 threshold fixed at 2.3.
legend(100, -200, legend=c("linear", "logarithmic", "exponential", "quadratic", "sine", "a-nonlinear", "BF01threshold=2.3"), 
       col = c("red", "darkblue", "magenta", "orange", "pink", "black", "green"), lty=1:2, cex=0.7)

##########################################################################################

#Effect Sizes (beta)

betas <- seq(-3, 3, length.out = 10) #beta sequence
nonlinear_beta <-c() #vector for the results
for (i in 1:10){
  nonlinear_beta <- cbind(nonlinear_beta, non_linear_test(beta = betas[i]))
}

beta_means <- colMeans(nonlinear_beta)
plot(seq(from =-3, to = 3, length = 10), beta_means[c(1, 7, 13, 19, 25, 31, 37, 43, 49, 55)], ylim = c(min(beta_means),max(beta_means)),
     type = "l", col = "red", main = "log(BF01) for different Beta values", 
     ylab = "log(BF01)", xlab = "Beta") #linear
lines(seq(from =-3, to = 3, length = 10), 
      beta_means[c(2, 8, 14, 20, 26, 32, 38, 44, 50, 56)], col = "darkblue") #logarithmic
lines(seq(from =-3, to = 3, length = 10), 
      beta_means[c(3, 9, 15, 21, 27, 33, 39, 45, 51, 57)], col = "magenta") #exponential
lines(seq(from =-3, to = 3, length = 10), 
      beta_means[c(4, 10, 16, 22, 28, 34, 40, 46, 52, 58)], col= "orange") #quadratic
lines(seq(from =-3, to = 3, length = 10), 
      beta_means[c(5, 11, 17, 23, 29, 35, 41, 47, 53, 59)], col = "pink") #sine
lines(seq(from =-3, to = 3, length = 10), 
      beta_means[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60)]) #a-nonlinear 
lines(seq(from =-3, to = 3, length = 10), y = c(rep(2.3, 10)),  col = "green", lwd = 2, lty = 2) #BF01 threshold fixed at 2.3.
legend(-0.75, -30, legend=c("linear", "logarithmic", "exponential", "quadratic", "sine", "a-nonlinear", "BF01threshold=2.3"), 
       col = c("red", "darkblue", "magenta", "orange", "pink", "black", "green"), lty=1:2, cex=0.7)



##########################################################################################

#Error Rate (sigma2)
sigma2 <- c(seq(from = 0.01,  to = 1, length.out = 10)) #sigma2
B0_sigma <- c() #vector for the results
for(i in 1: 10){
  sigmax <- non_linear_test(sigma2 = sigma2[i])
  B0_sigma <- cbind(B0_sigma, sigmax)
}

sigma2_means <- colMeans(B0_sigma)
plot(sigma2, y=sigma2_means[c(1, 7, 13, 19, 25, 31, 37, 43, 49, 55)], ylim = c(min(sigma2_means), max(sigma2_means)), 
     type = "l", col = "red",main ="log(BF01) for different sigma values", 
     ylab = "log(BF01)", xlab = "sigma2") #linear
lines(sigma2, y = sigma2_means[c(2, 8, 14, 20, 26, 32, 38, 44, 50, 56)], col = "darkblue") #logarithmic
lines(sigma2, y = sigma2_means[c(3, 9, 15, 21, 27, 33, 39, 45, 51, 57)], col = "magenta") #exponential
lines(sigma2, y = sigma2_means[c(4, 10, 16, 22, 28, 34, 40, 46, 52, 58)], col = "orange") #quadratic 
lines(sigma2, y= sigma2_means[c(5, 11, 17, 23, 29, 35, 41, 47, 53, 59)], col = "pink") #sine 
lines(sigma2, y = sigma2_means[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60)]) #a-nonlinear
lines(sigma2, y = c(rep(2.3, 10)), col = "green", lwd = 2, lty = 2) #BF01 threshold fixed at 2.3. 
legend(0.5, -40,
       legend=c("linear", "logarithmic", "exponential", "quadratic", "sine", "a-nonlinear", "BF01threshold=2.3"), 
       col = c("red", "darkblue", "magenta", "orange", "pink", "black", "green"), lty=1:2, cex=0.7)

##########################################################################################


##########################################################################################
                                    #COMPLEX MODEL TEST
##########################################################################################


#Matrix for average Bayes Factors for each model against all other models 
mean_BFs <- matrix(NA, nrow = 500, ncol = 27)
colnames(mean_BFs) <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", 
                        "M10","M11", "M12", "M13","M14", "M15", "M16", "M17", "M18",
                        "M19", "M20", "M21", "M22", "M23", "M24", "M25", "M26", "M27")
complex_test <- function(n = 100, sigma2 = 0.1, k=4, b=2){#default n=100 and sigma2=0.1
  for (j in 1:500){
    
    features <- cbind("x1" = rnorm(n, 0, .5), "x2" = rnorm(n, 0, .5), "x3" = rnorm(n, 0, .5))
    beta <- c(b, b, 0) #discuss how to adjust beta values 
    error <- rnorm(n,sd=sqrt(sigma2))
    y <- beta[1]*exp(features[,1]) + beta[2]*features[,2] + beta[3]*features[,3] + error
    data <- as.data.frame(cbind(y, features))
    
    #model fits & BIC calculations
    bics <- c(rep(NA, 27)) #vector for bic scores of the models fitted below
    
    M1 <- gam(y ~ 1,data=data); bics[1]=(-2)*head(logLik(M1))+2*log(n)
    
    M2 <- gam(y ~ 1 + x1, data=data); bics[2]=(-2)*head(logLik(M2))+3*log(n)
    M3 <- gam(y ~ 1 + x2, data=data); bics[3]=(-2)*head(logLik(M3))+3*log(n)
    M4 <- gam(y ~ 1 + x3, data=data); bics[4]=(-2)*head(logLik(M4))+3*log(n)
    
    M5 <- gam(y ~ 1 + s(x1, k=k), data=data);bics[5]=(-2)*head(logLik(M5))+5*log(n)
    M6 <- gam(y ~ 1 + s(x2, k=k), data=data);bics[6]=(-2)*head(logLik(M6))+5*log(n) 
    M7 <- gam(y ~ 1 + s(x3, k=k), data=data);bics[7]=(-2)*head(logLik(M7))+5*log(n)
    
    M8 <- gam(y ~ 1 + x1 + x2, data=data); bics[8]=(-2)*head(logLik(M8))+4*log(n)
    M9 <- gam (y ~ 1 + x1 + x3, data=data); bics[9]=(-2)*head(logLik(M9))+4*log(n)
    M10 <- gam(y ~ 1 + x2 + x3, data=data); bics[10]=(-2)*head(logLik(M10))+4*log(n)
    
    M11 <- gam(y ~ 1 + s(x1, k=k) + x2, data=data); bics[11]=(-2)*head(logLik(M11))+6*log(n) #true
    M12 <- gam(y ~ 1 + s(x1, k=k) + x3, data=data); bics[12]=(-2)*head(logLik(M12))+6*log(n)
    
    M13 <- gam(y ~ 1 + x1 + s(x2, k=k), data=data); bics[13]=(-2)*head(logLik(M13))+6*log(n)
    M14 <- gam(y ~ 1 + x1 + s(x3, k=k), data=data); bics[14]=(-2)*head(logLik(M14))+6*log(n)
    
    M15 <- gam(y ~ 1 + s(x2, k=k) + x3, data=data); bics[15]=(-2)*head(logLik(M15))+6*log(n)
    M16 <- gam(y ~ 1 + x2 + s(x3, k=k), data=data); bics[16]=(-2)*head(logLik(M16))+6*log(n)
    
    M17 <- gam(y ~ 1 + s(x1, k=k) + s(x2, k=k), data=data); bics[17]=(-2)*head(logLik(M17))+8*log(n)
    M18 <- gam(y ~ 1 + s(x1, k=k) + s(x3, k=k), data=data); bics[18]=(-2)*head(logLik(M18))+8*log(n)
    M19 <- gam(y ~ 1 + s(x2, k=k) + s(x3, k=k), data=data); bics[19]=(-2)*head(logLik(M19))+8*log(n)
    
    M20 <- gam(y ~ 1 + x1 + x2 + x3, data = data); bics[20]=(-2)*head(logLik(M20))+5*log(n)
    
    M21 <- gam(y ~ 1 + s(x1, k=k) + x2 + x3, data=data); bics[21]=(-2)*head(logLik(M21))+7*log(n)
    M22 <- gam(y ~ 1 + x1 + s(x2, k=k) + x3, data=data); bics[22]=(-2)*head(logLik(M22))+7*log(n)
    M23 <- gam(y ~ 1 + x1 + x2 + s(x3, k=k), data=data); bics[23]=(-2)*head(logLik(M23))+7*log(n)
    
    M24 <- gam(y ~ 1 + s(x1, k=k) + s(x2, k=k) + x3, data=data); bics[24]=(-2)*head(logLik(M24))+9*log(n)
    M25 <- gam(y ~ 1 + s(x1, k=k) + x2 + s(x3, k=k), data=data); bics[25]=(-2)*head(logLik(M25))+9*log(n)
    M26 <- gam(y ~ 1 + x1 + s(x2, k=k) + s(x3, k=k), data=data); bics[26]=(-2)*head(logLik(M26))+9*log(n)
    
    M27 <- gam(y ~ 1 + s(x1, k=k) + s(x2, k=k) + s(x3, k=k), data=data); bics[27]=(-2)*head(logLik(M27))+11*log(n)
    
    #matrix for Bayes Factors 
    BFs <- matrix(NA, nrow = 27, ncol = 27)
    colnames(BFs) <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", 
                       "M9", "M10","M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18",
                       "M19", "M20", "M21", "M22", "M23", "M24", "M25", "M26", "M27")
    rownames(BFs) <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", 
                       "M9", "M10","M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18",
                       "M19", "M20", "M21", "M22", "M23", "M24", "M25", "M26", "M27")
    
    
    for (k in 1:27) {
      bf_vec<-c()
      for (x in 1:26) {
        bf <- log(exp((bics[-k][x] - bics[k]) /2)) 
        bf_vec <- c(bf_vec, bf)
      }
      bf_vec <- R.utils::insert(bf_vec, ats=k, values=0) #inserting 0 to the diagonals
      BFs[, k] <- round(bf_vec, digits=3)
    }
    
    BFs[BFs == "-Inf"]  <- -1e5; BFs[BFs == "Inf"]  <- 1e5 #to deal with infinite values
    mean_BFs[j,] <- colMeans(BFs)
    
    result <- list("the best model chosen" = sort(colMeans(mean_BFs))[27],
                   "Mean Bayes factors for the models (unsorted" = colMeans(mean_BFs),
                   "Mean Bayes factors for the models (sorted)" = sort(colMeans(mean_BFs)),
                   "Bayes Factors for each trial" = mean_BFs)
  }
  return(result)
}

##########################################################################################

#Sample Sizes (n)

n_vec_BF <- c(rep(NA, 4)) #for the BF value of the best model chosen
n_vec_name <- c(rep(NA, 4)) #exact name of the best model chosen
n_sizes <- c(100, 250, 500, 1000)
for (i in 1:4){
  x <- complex_test(n = n_sizes[i])
  n_vec_BF[i] <- x[[1]] 
  n_vec_name[i] <- names(x[[1]])
}
names(n_vec_BF) <- n_vec_name

##########################################################################################

#Effect Sizes (beta)
betas <- seq(-3, 3, length.out = 10) #beta sequence
beta_vec_BF <- c(rep(NA, 10)) #for the BF value of the best model chosen
beta_vec_name <- c(rep(NA, 10)) #exact name of the best model chosen

for (i in 1:10){
  x<-complex_test(b = betas[i])
  beta_vec_BF[i] <- x[[1]]
  beta_vec_name[i] <- names(x[[1]])
}

names(beta_vec_BF) <- beta_vec_name

##########################################################################################
##########################################################################################
##########################################################################################





