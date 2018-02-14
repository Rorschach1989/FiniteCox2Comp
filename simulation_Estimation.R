#######################################################
# A Brief Illustration of Estimating a 2-component
# Mixture-of-Experts Proportional Hazard Regression
# Model
#######################################################
#######################################################
source('EMEstimation.R')

#######################################################
# Specify some global parameters, like number of sample
# points, dimension of Response Covariates and dimension
# of Subgroup related Covariates
N <- 500
p <- 2
q <- 2
#######################################################
# Since this is a simulation study, we specify the true
# parameters here
beta10 <- c(1,1)
beta20 <- c(1.5,1.5)
gamma0 <- c(0.5,-1.5)

M <- 1000
count <- 0
result_Est <- matrix(0,M,(2*p+q))
result_Var <- matrix(0,M,(2*p+q))
nonPSD_failure <- 0

set.seed(7777777)
while (TRUE){
  # Covariates
  z1 <- rnorm(N,mean=1,sd=1)
  z2 <- rbinom(N,1,0.5)
  Z <- cbind(z1 = z1,z2 = z2)
  x1 <- runif(N,-3,3)
  X <- cbind(x0 = 1,x1 = x1)
  SRate <- exp(X%*%gamma0)/(1+exp(X%*%gamma0))
  # Generate Data
  ksi <- rbinom(N,1,SRate)
  rate1 <- exp(Z%*%(beta10 + beta20))
  rate2 <- exp(Z%*%beta10)
  pop1 <- rexp(N,rate=rate1)
  pop2 <- rexp(N,rate=rate2)
  pop <- ksi*pop1+(1-ksi)*pop2
  # This censoring mechansim controls the censoring rate around 10-15%
  censor <- rexp(N,rate=rate2) + .3
  y <- pmin(pop,censor)
  delta <- as.numeric(pop <= censor)
  DataA <- data.frame(y,delta,Z,X)
  if (count >= M) break
  
  # config for which part of the dataset would be used for analysis
  config <- list(y = 'y', delta = 'delta', Z = c("z1", "z2"), X = c("x0", "x1"))
  # The variance estimation of the estimator is performed either in a 
  # bootstrap fashion, or using the full semiparametric Fisher Information (in some loose sence)
  # bootstrap is extremely slow while more robust against small sample size and (perhaps?) some
  # mysterious level of model misspecification, if speed is of concern, set bootR = 0 in the 
  # estimation step
  
  # Here true value is chosen as initial value which is unrealistic, this is for 
  # preventing osciliations between multiple global maxmimums (ME style Finite Mixture is 
  # identified only up to transformations if not constrained in a certain way, while making
  # the constrainits explicit during optimization is not appropriate since in real world 
  # problems finding the best combination for interpretation would be easy
  betaInitVal <- matrix(c(beta10,beta20),1,2*p)
  gammaInitVal <- matrix(gamma0,1,q)
  estResult <- mixEst(DataA,config,
                      betaInitVal,gammaInitVal,precision = 1e-7,bootR = 0)
  if (!(inherits(estResult$SSE,"error") | inherits(estResult$SSE,"warning"))){
    result_Est[count,] <- c(estResult$beta_Est, estResult$gamma_Est)
    #fill in the variance estimation result matrix
    result_Var[count,] <- estResult$SSE
    count <- count + 1
    cat('\r','Progress Counter: ', toString(count),'/',toString(M))
  }else{
    nonPSD_failure <- nonPSD_failure + 1
    next
  }
}

cat('\n')
print('==================Evaluations==================')
print('Bias:')
print(apply(result_Est, 2, mean) - c(beta10, beta20, gamma0))
print('SE:')
print(apply(result_Est, 2, sd))
print('SEE')
print(apply(result_Var, 2, mean))





