source('EMbasic.R')
source('FisherInfoMat.R')
library(boot)

mixEst <- function(DataR,config,
                   betaInit,gammaInit,precision = 1e-5,bootR = 499){
  # Parse dataframe w.r.t. config
  y <- config$y
  # reorder input dataframe w.r.t. min(event_time, censor_time)
  DataR <- DataR[order(DataR[[y]]),]
  Z <- as.matrix(subset(DataR, select = config$Z))
  X <- as.matrix(subset(DataR, select = config$X))
  delta <- as.vector(DataR[[config$delta]])
  p <- ncol(Z)
  q <- ncol(X)
  #Estimation Phase
  #betaInit and gammaInit should be initial value matrices with same row length, each row represent an initial value choice
  Runs <- nrow(betaInit)
  likVec <- list()
  iMax <- 1
  for (i in 1:Runs){
    rateTmp <- exp(X%*%gammaInit[i,])/(1+exp(X%*%gammaInit[i,]))
    likVec[[i]] <- EMfull(Z,X,delta,
                          betaInit[i,],profR(betaInit[i,],Z,delta,rateTmp),gammaInit[i,],precision = precision)
    if (i == 1){
      maxlik = likVec[[i]]$likConv
    }else{
      if (maxlik < likVec[[i]]$likConv){
        maxlik <- likVec[[i]]$likConv
        iMax <- i
      }
    }
  }
  estRaw <- likVec[[iMax]]
  rm(likVec)
  etaFinal <- estepCalc(estRaw$beta, estRaw$hazard, estRaw$gamma, Z, X, delta)
  # Variance estimation phase, bootstrap or asymptotic
  if (bootR){
    forBoot <- function(DataR_){
      # Parse dataframe w.r.t. config
      y <- config$y
      # reorder input dataframe w.r.t. min(event_time, censor_time)
      DataR_ <- DataR_[order(DataR[[y]]),]
      Z <- as.matrix(subset(DataR_, select = config$Z))
      X <- as.matrix(subset(DataR_, select = config$X))
      delta <- as.vector(DataR_[[config$delta]])
      p <- ncol(Z)
      q <- ncol(X)
      EMfull(Z,X,delta,
             estRaw$beta,HaZ = 1,gamma1 = estRaw$gamma, precision = precision, boot = TRUE)
    }
    bootResult <- censboot(DataR,forBoot,R = bootR)
    varRaw <- bootResult$t
    for (j in 1:q){
      varRaw <- varRaw[which(abs(varRaw[,(2*p+j)]) < 36),]
    }
    ssE <- apply(varRaw,2,sd)
    return(list(beta_Est = estRaw$beta, gamma_Est = estRaw$gamma, haz_Est = estRaw$hazard, SSE = ssE))
  }else{
    beta1Est <- estRaw$beta[1:p]
    beta2Est <- estRaw$beta[(p+1):(2*p)]
    gammaEst <- estRaw$gamma
    haz_Est = estRaw$hazard
    AsymtHess <- FisherInfo(beta1 = beta1Est, beta2 = beta2Est, gamma1 = gammaEst, 
                            hazard = haz_Est, 
                            eta = as.vector(etaFinal), Z = Z, X = X)
    varEstTrial <- tryCatch(hechol <- solve(chol(AsymtHess)),
                            error = function(e) e,
                            warning = function(w) w)
    if ((inherits(varEstTrial,"error") | inherits(varEstTrial,"warning"))){
      return(list(beta_Est = c(beta1Est,beta2Est), gamma_Est = gammaEst, haz_Est = haz_Est, SSE = varEstTrial))
    }else{
      sdEst <- sqrt(diag(hechol %*% t(hechol))[1:(2*p + q)])
      return(list(beta_Est = c(beta1Est,beta2Est), gamma_Est = gammaEst, haz_Est = haz_Est, SSE = sdEst))
    }
  }
}
