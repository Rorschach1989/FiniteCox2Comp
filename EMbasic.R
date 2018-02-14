source('utils.R')
#Z should be in a pre-sorted order
mstepforBeta <- function(betaT,Z,delta,Eta){
  beta1 <- betaT[1:p]
  beta2 <- betaT[(p+1):(2*p)]
  bZ1 <- Z%*%(beta1 + beta2)
  bZ2 <- Z%*%beta1
  ebZ1 <- exp(bZ1)
  ebZ2 <- exp(bZ2)
  -sum(delta*(Eta*bZ1 + (1-Eta)*bZ2 - log(rcumsumr(Eta*ebZ1 + (1-Eta)*ebZ2))))
}

gradmstepforBeta <- function(betaT,Z,delta,Eta){
  beta1 <- betaT[1:p]
  beta2 <- betaT[(p+1):(2*p)]
  bZ1 <- Z%*%(beta1 + beta2)
  bZ2 <- Z%*%beta1
  ebZ1 <- exp(bZ1)
  ebZ2 <- exp(bZ2)
  deNom <- as.vector(rcumsumr(Eta*ebZ1 + (1-Eta)*ebZ2))
  nom1 <- apply(as.vector(Eta*ebZ1)*Z,2,rcumsumr)
  nom2 <- apply(as.vector((1-Eta)*ebZ2)*Z,2,rcumsumr)
  grad1 <- apply(delta*Z,2,sum) - apply(delta/deNom*(nom1 + nom2),2,sum)
  grad2 <- apply(delta*Eta*Z,2,sum) - apply(delta/deNom*nom1,2,sum)
  -c(grad1,grad2)
}

profR <- function(betaT,Z,delta,Eta){
  beta1 <- betaT[1:p]
  beta2 <- betaT[(p+1):(2*p)]
  bZ1 <- Z%*%(beta1 + beta2)
  bZ2 <- Z%*%beta1
  ebZ1 <- exp(bZ1)
  ebZ2 <- exp(bZ2)
  deNom <- as.vector(rcumsumr(Eta*ebZ1 + (1-Eta)*ebZ2))
  delta/deNom
}

mstepforGamma <- function(gamma1,Eta,X){
  tmp <- X%*%gamma1
  -sum(Eta*tmp - log(1+exp(tmp)))
}

gradmstepforGamma <- function(gamma1,Eta,X){
  tmp <- exp(X%*%gamma1)
  -apply(as.vector(Eta-tmp/(1+tmp))*X,2,sum)
}



phiSeries <- function(betaR,HaZ,Z,delta){
  bZ <- Z%*%betaR
  ebZ <- exp(bZ)
  cH <- cumsum(HaZ)
  (HaZ*ebZ)^(delta)*exp(-cH*ebZ)
}

estepCalc <- function(betaT,HaZ,gamma1,Z,X,delta){
  beta1 <- betaT[1:p]
  beta2 <- betaT[(p+1):(2*p)]
  mixprob <- exp(X%*%gamma1)/(1+exp(X%*%gamma1))
  s1 <- phiSeries(beta1+beta2,HaZ,Z,delta)*mixprob
  s2 <- phiSeries(beta1,HaZ,Z,delta)*(1-mixprob)
  s1/(s1+s2)
}

loglikFull <- function(betaT,HaZ,gamma1,Z,X,delta){
  beta1 <- betaT[1:p]
  beta2 <- betaT[(p+1):(2*p)]
  mixprob <- exp(X%*%gamma1)/(1+exp(X%*%gamma1))
  s1 <- phiSeries(beta1+beta2,HaZ,Z,delta)*mixprob
  s2 <- phiSeries(beta1,HaZ,Z,delta)*(1-mixprob)
  sum(log(s1+s2))
}

EMfull <- function(Z,X,delta,
                   betaT,HaZ,gamma1,precision = 1e-5,
                   boot = FALSE,eta = 0,maxiter = 0){
  # full EM algorithm in a single shot
  # eta is only needed for bootstrapping
  err <- 5
  i <- 0
  if (boot) {
    HaZ <- profR(betaT,Z,delta,eta)
  }
  loglik <- loglikFull(betaT,HaZ,gamma1,Z,X,delta)
  while (TRUE){
    i <- i+1
    #E-step
    eta <- as.vector(estepCalc(betaT,HaZ,gamma1,Z,X,delta))
    if (sum(is.nan(eta)) > 0){
      # handler for numerical failure encountered at initialization phase, rather crude
      eta <- rep(.5,N)
    }
    vartmp <- c(betaT,gamma1)
    #M-step, for beta
    betaT <- optim(betaT,mstepforBeta,gradmstepforBeta,method = "BFGS",
                   Eta = eta,Z = Z,delta = delta)$par
    HaZ <- profR(betaT,Z,delta,eta)
    gamma1 <- optim(gamma1,mstepforGamma,gradmstepforGamma,method = "BFGS",
                    Eta = eta,X = X)$par
    tmplik <- loglik
    loglik <- loglikFull(betaT,HaZ,gamma1,Z,X,delta)
    err <- lpNorm(vartmp-c(betaT,gamma1),2)
    if (loglik - tmplik < precision & err < precision) break
    if (maxiter & maxiter < i) break
  }
  if (boot) {
    return(c(betaT,gamma1))
  }else{
    return(list(beta = betaT,hazard = HaZ,gamma = gamma1,steps = i,likConv = loglik, tolerance = err, Eta = eta))
  }
}

EMTestBasic <- function(Z,X,delta,
                        betaT,HaZ,gamma1,precision = 1e-5,maxiter = 0){
  # partial EM algorithm used for EM-test
  err <- 5
  i <- 0
  loglik <- loglikFull(betaT,HaZ,gamma1,Z,X,delta)
  while (TRUE){
    i <- i+1
    #E-step
    eta <- as.vector(estepCalc(betaT,HaZ,gamma1,Z,X,delta))
    if (sum(is.nan(eta)) > 0){
      eta <- rep(.5,N)
    }
    vartmp <- betaT
    #M-step, for beta
    betaT <- optim(betaT,mstepforBeta,gradmstepforBeta,method = "BFGS",
                   Eta = eta,Z = Z,delta = delta)$par
    HaZ <- profR(betaT,Z,delta,eta)
    tmplik <- loglik
    loglik <- loglikFull(betaT,HaZ,gamma1,Z,X,delta)
    err <- lpNorm(vartmp-betaT,2)
    if (loglik - tmplik < precision & err < precision) break
    if (maxiter & maxiter < i) break
  }
  list(beta = betaT, hazard = HaZ)
}