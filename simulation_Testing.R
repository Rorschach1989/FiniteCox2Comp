#######################################################
# A Brief Illustration of EM-Test for checking 
# heterogeneity in a potentially 2 component Cox Model
#######################################################
#######################################################
source('EMTest.R')

#######################################################
# Specify some global parameters, like number of sample
# points, dimension of Response Covariates and dimension
# of Subgroup related Covariates
N <- 100
p <- 2
q <- 2
#######################################################
# Since this is a simulation study, we specify the true
# parameters here, the ground truth is set to be 
# homogeneous, which requires the proposed algorithm
# to control type I error
beta10 <- c(1,1)
beta20 <- c(0,0)
gamma0 <- c(0.5,-1.5)

M <- 100
count <- 0

# Boostrap configuration, resample 99 times each
bootR <- 99
eSeq <- c(0,3,9)
result_Test_Boot <- matrix(0, M, length(eSeq)) 

# EM-Test configuration
J <- 2
gammaJ <- matrix(runif(J*q),J,q)
print('gammaJ configuration: ')
print(gammaJ)

start <- Sys.time()
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
  testResult <- EMTest(DataR = DataA, config = config, gammaJ = gammaJ, eSeq = eSeq, bootR = bootR, precision = 1e-7)
  result_Test_Boot[count,] <- testResult$Pvalue_boot
  count <- count + 1
  cat('\r','Progress Counter: ', toString(count),'/',toString(M), 'Time elapsed', Sys.time() - start)
}

cat('\n')
print('==================Evaluations==================')
print('Type-I Error @ 0.05 significance:')
print(apply(result_Test_Boot, 2, function(x) sum(x < .05)/M))
print('Type-I Error @ 0.10 significance:')
print(apply(result_Test_Boot, 2, function(x) sum(x < .10)/M))


