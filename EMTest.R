source('EMbasic.R')
library(survival)
library(boot)

EMTestSinle <- function(Z,X,delta,
                        betaT,HaZ,loglik0,
                        gamma1,precision = 1e-5,eSeq){
  EM_0 <- EMTestBasic(Z,X,delta,
                      betaT,HaZ,
                      gamma1,precision = precision)
  loglik1 <- loglikFull(EM_0$beta,EM_0$hazard,gamma1,Z,X,delta)
  #create a vector to store the result
  l <- length(eSeq)
  testResult <- rep(0,l)
  nameVec <- rep("",l)
  num_last <- 0
  for (i in 1:l){
    if (eSeq[i] == 0){
      testResult[i] <- 2*(loglik1 - loglik0)
      nameVec[i] <- paste0("EM_",toString(eSeq[i]))
      next
    }
    for (j in 1:(eSeq[i] - num_last)){
      #perform several "complete" EM-steps
      eta <- as.vector(estepCalc(betaT,HaZ,gamma1,Z,X,delta))
      betaT <- optim(betaT,mstepforBeta,gradmstepforBeta,method = "BFGS",
                     Eta = eta,Z = Z,delta = delta)$par
      HaZ <- profR(betaT,Z,delta,eta)
      gamma1 <- optim(gamma1,mstepforGamma,gradmstepforGamma,method = "BFGS",
                      Eta = eta,X = X)$par
    }
    EM_1 <- EMTestBasic(Z,X,delta,
                        betaT,HaZ,gamma1,precision = precision)
    loglik1 <- loglikFull(EM_1$beta,EM_1$hazard,gamma1,Z,X,delta)
    testResult[i] <- 2*(loglik1 - loglik0)
    num_last <- eSeq[i]
    nameVec[i] <- paste0("EM_",toString(eSeq[i]))
  }
  names(testResult) <- nameVec
  testResult
}

EMTestAsym <- function(DataR,config,
                   gammaJ,precision = 1e-5,eSeq){
  # Method for calculating EM-Test statistic
  # Parse dataframe w.r.t. config
  y <- config$y
  # reorder input dataframe w.r.t. min(event_time, censor_time)
  DataR <- DataR[order(DataR[[y]]),]
  Z <- as.matrix(subset(DataR, select = config$Z))
  X <- as.matrix(subset(DataR, select = config$X))
  delta <- as.vector(DataR[[config$delta]])
  p <- ncol(Z)
  q <- ncol(X)
  # Perform ordinary Cox model assuming no subgroup effect
  DataSubset <- subset(DataR, select = c(config$y, config$delta, config$Z))
  formula_ <- paste0('Surv(', config$y, ",", config$delta, ')~.')
  coxEst <- coxph(as.formula(formula_),data = DataSubset)
  betaT <- c(coxEst$coefficients,rep(0,p))
  HaZ <- profR(betaT,Z,delta,Eta = 1)
  #EM(0) should be always calculated, but may be omitted in the final result if eSeq does not 
  #contain 0
  if (!is.matrix(gammaJ)){
    stop('gammaJ should be a matrix with rows containing consecutive gamma tries')
  }
  if (!ncol(gammaJ) == q){
    stop('gamma tries should match the dimension of X')
  }
  fullResult <- matrix(0, nrow(gammaJ), length(eSeq))
  for (j in 1:nrow(gammaJ)){
    gamma_ = gammaJ[j,]
    loglik0 <- loglikFull(betaT,HaZ,gamma_,Z,X,delta)
    fullResult[j,] <- EMTestSinle(
      Z,X,delta,
      betaT,HaZ,loglik0,
      gamma1 = gamma_, precision = precision, eSeq = eSeq
    )
  }
  apply(fullResult, 2, max)
}

EMTestBoot <- function(DataR,config,
                       gammaJ,eSeq,bootR,precision = 1e-5){
  # Parse dataframe w.r.t. config
  y <- config$y
  # reorder input dataframe w.r.t. min(event_time, censor_time)
  DataR <- DataR[order(DataR[[y]]),]
  Z <- as.matrix(subset(DataR, select = config$Z))
  X <- as.matrix(subset(DataR, select = config$X))
  delta <- as.vector(DataR[[config$delta]])
  p <- ncol(Z)
  q <- ncol(X)
  # Perform ordinary Cox model assuming no subgroup effect
  DataSubset <- subset(DataR, select = c(config$y, config$delta, config$Z))
  formula_ <- paste0('Surv(', config$y, ",", config$delta, ')~.')
  coxEst <- coxph(as.formula(formula_),data = DataSubset)
  coxSurv <- survfit(coxEst)
  DataSubset_ <- DataSubset
  DataSubset_[[config$delta]] <- 1 - DataSubset_[[config$delta]]
  coxEstCens <- coxph(as.formula(formula_),data = DataSubset_)
  coxSurvCens <- survfit(coxEstCens)
  coxBoot <- censboot(DataR, EMTestAsym, R=bootR,
                      F.surv = coxSurv,
                      G.surv = coxSurvCens,
                      cox = coxEst,
                      sim = 'model',
                      config = config, gammaJ = gammaJ, eSeq = eSeq
                      )
  coxBoot
}

EMTest <- function(DataR, config,
                   gammaJ, eSeq,
                   bootR, precision = 1e-5){
  chisq_df <- length(config$X)
  if (nrow(gammaJ) == 1){
    if (bootR){
      testResult <- EMTestBoot(DataR,config,
                               gammaJ,eSeq,bootR,precision)
      Test_Statistic <- testResult$t0
      Pvalue_asym <- NULL
      Pvalue_boot <- apply(testResult$t - testResult$t0, 2, function(x) sum(x >= 0)/(bootR + 1))
    }else{
      testResult <- EMTestAsym(DataR,config,
                               gammaJ,eSeq,bootR,precision)
      Test_Statistic <- testResult
      Pvalue_asym <- pchisq(Test_Statistic, df = chisq_df)
      Pvalue_boot <- NULL
    }
  }else{
    if(!bootR){
      stop('In a multi-sequence EM-Test, non-boostrap method are unsupported')
    }
    testResult <- EMTestBoot(DataR,config,
                             gammaJ,eSeq,bootR,precision)
    Test_Statistic <- testResult$t0
    Pvalue_asym <- NULL
    Pvalue_boot <- apply(testResult$t - testResult$t0, 2, function(x) sum(x >= 0)/(bootR + 1))
  }
  return(list(
    Test_Statistic = Test_Statistic, Pvalue_asym = Pvalue_asym, Pvalue_boot = Pvalue_boot
  ))
}