source('utils.R')
FisherInfo <- function(beta1,beta2,gamma1,hazard,
                       eta,Z,X){
  #seperate parameters
  index <- which(hazard * 1e1000 != 0)
  h <- hazard[index]
  ne <- length(h)
  H <- cumsum(hazard)
  ebZ1 <- as.vector(exp(Z%*%(beta1 + beta2)))
  ebZ2 <- as.vector(exp(Z%*%beta1))
  ebZH1 <- as.vector(ebZ1*H)
  ebZH2 <- as.vector(ebZ2*H)
  Xg <- exp(X%*%gamma1)
  #STEP I E(SS^T|obs)
  #coefficient construction phase
  coefb1 <- -ebZH1
  coefb2 <- -ebZH2
  coefh1 <- ebZ1
  coefh2 <- ebZ2
  coefg <- as.vector(Xg/(1+Xg))
  #index correction, very important
  coefb1[index] <- coefb1[index] + 1
  coefb2[index] <- coefb2[index] + 1
  #coefh1[-index] <- 0
  #coefh2[-index] <- 0
  #vector filling phase
  lb1 <- coefb1*Z
  lb2 <- coefb2*Z
  lg <- X
  #block construction phase
  tmpb1 <- t(sqrt(eta)*lb1)%*%(sqrt(eta)*lb1) - t(eta*lb1)%*%(eta*lb1)
  tmpb2 <- t(sqrt(1-eta)*lb2)%*%(sqrt(1-eta)*lb2) - t((1-eta)*lb2)%*%((1-eta)*lb2)
  tmpb3 <- -t((1-eta)*lb2)%*%(eta*lb1)
  h1b1b1 <- tmpb1 + tmpb2 + tmpb3 + t(tmpb3)
  h1b1b2 <- tmpb1 + tmpb3
  h1b2b2 <- tmpb1
  rm(tmpb1,tmpb2,tmpb3)
  tmpb1 <- t(sqrt(eta)*lb1)%*%(sqrt(eta)*lg) - t(eta*lb1)%*%(eta*lg)
  tmpb2 <- -t((1-eta)*lb2)%*%(eta*lg)
  tmpg <- t(sqrt(eta)*lg)%*%(sqrt(eta)*lg) - t(eta*lg)%*%(eta*lg)
  h1gg <- tmpg
  h1b1g <- tmpb1 + tmpb2 # Cb1%*%t(Cgamma) should be zero
  h1b2g <- tmpb1 # Cb2%*%t(Cgamma) should be zero
  rm(tmpb1,tmpb2,tmpg)
  #h part is crucial
  s1 <- -(eta-eta^2)*(coefh1 - coefh2)*lb1
  s2 <- (eta-eta^2)*(coefh1 - coefh2)*lb2
  s3 <- -(eta-eta^2)*(coefh1 - coefh2)*lg
  tmpb1 <- apply(s1,2,rcumsumr)[index,]
  tmpb2 <- apply(s2,2,rcumsumr)[index,]
  tmpg <- apply(s3,2,rcumsumr)[index,]
  h1b1h <- (tmpb1 + tmpb2)#[index,] #(Ch%*%t(Cb1)) is zero
  h1b2h <- (tmpb1)#[index,] #(Ch%*%t(Cb2)) is zero
  h1gh <- (tmpg)#[index,] #(Ch + Chstar)%*%t(Cgamma) is zero
  element <- rcumsumr((eta-eta^2)*(coefh1 - coefh2)^2)[index]
  h1hh <- sapply(1:ne,nondiagpick,x=element)
  h1hh <- h1hh + t(h1hh)
  #assemble matrices
  hess1 <- rbind(cbind(h1b1b1,h1b1b2,h1b1g,t(h1b1h)),cbind(t(h1b1b2),h1b2b2,h1b2g,t(h1b2h)),
                 cbind(t(h1b1g),t(h1b2g),h1gg,t(h1gh)),cbind(h1b1h,h1b2h,h1gh,h1hh))
  rm(tmpb1,tmpb2,tmpg,h1b1b1,h1b1b2,h1b2b2,h1b1g,h1b2g,h1gg,h1b1h,h1b2h,h1gh,h1hh)
  #STEP II E(H|obs)
  #beta_2 block is just part of the beta_1 block, so calculate it first
  #beta_2 - beta_2
  b2vec <- ebZH1
  Ztil2 <- sqrt(b2vec*eta)*Z
  h2b2b2 <- -t(Ztil2)%*%Ztil2
  #beta_2 - beta_1 is just the same!
  h2b1b2 <- h2b2b2
  #beta_1 - beta_1
  b1vec <- ebZH2
  Ztil1 <- sqrt(b1vec*(1-eta))*Z
  h2b1b1 <- -t(Ztil1)%*%Ztil1 + h2b2b2
  #beta_1-gamma,beta_2-gamma are 0s
  h2b1g <- h2b2g <- matrix(0,q,p)
  #gamma-gamma
  gvec <- as.vector(Xg/((1+Xg)^2))
  Xtil <- sqrt(gvec)*X
  h2gg <- -t(Xtil)%*%Xtil
  #gamma-h are 0s
  h2gh <- matrix(0,ne,q)
  #beta_2-h block is also part of the beta_1-h block
  h2b2h <- -apply(eta*coefh1*Z,2,rcumsumr)[index,]
  #beta_1-h
  h2b1h <- -apply((1-eta)*coefh2*Z,2,rcumsumr)[index,] + h2b2h
  #h-h
  h2hh <- -diag(1/(h^2))
  hess2 <- rbind(cbind(h2b1b1,t(h2b1b2),t(h2b1g),t(h2b1h)),cbind(h2b1b2,h2b2b2,t(h2b2g),t(h2b2h)),
                 cbind(h2b1g,h2b2g,h2gg,t(h2gh)),cbind(h2b1h,h2b2h,h2gh,h2hh))
  rm(h2b1b1,h2b1b2,h2b2b2,h2b1g,h2b2g,h2gg,h2b1h,h2b2h,h2gh,h2hh)
  - hess2 - hess1
}