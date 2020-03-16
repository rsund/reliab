###############################################################################
#
# R-implementation for the calculation of Tarkkonen's rho (reli)
#   Tarkkonen L, Vehkalahti K (2005). Measurement errors in multivariate measurement scales.
#   Journal of Multivariate Analysis 96 (2005) 172–189, doi:10.1016/j.jmva.2004.09.007
#
# Implementations for calculating Alpha and Omega (reli)
#
# Functions (rlbsimul, simuexp) intended to conduct similar simulation studies
#   as reported in Vehkalahti (2000). http://hdl.handle.net/10138/10570
#
# Author(s):
# Reijo Sund (ORCID: 0000-0002-6268-8117), https://connect.uef.fi/en/reijo.sund
#
###############################################################################

reli <- function(R,B,W=B,O=diag(ncol(B)),M=rep(1,nrow(B)),sum=FALSE,scores=TRUE) {
  rhoE3 <- diag(diag(diag(t(W) %*% B %*% O %*% t(B) %*% W),ncol=ncol(B)) %*% solve(diag(diag(t(W) %*% R %*% W),ncol=ncol(B))))
  rhoE2 <- diag(solve(diag(ncol(B)) + diag(diag(t(W) %*% diag(diag(R-(B %*% t(B)))) %*% W),ncol=ncol(B)) %*% solve(diag(diag(t(W) %*% B %*% O %*% t(B) %*% W),ncol=ncol(B))) ))
  reli <- cbind(rhoE2,rhoE3)
  rownames(reli) <- colnames(W)

  RU <- R-diag(diag(R-(B %*% t(B))))
  oa <- diag(t(B) %*% RU %*% B)
  ob <- diag(t(B) %*% R %*% B)
  
  D <- diag(nrow=nrow(B))  # diag(M)
  S <- D %*% R %*% D
  dS <- diag(diag(S))
  p <- nrow(W)
  aa <- diag(t(W) %*% S %*% W)
  ab <- diag(t(W) %*% dS %*% W)
  
if (scores) {  
  rhoScore <- diag(t(B) %*% solve(R) %*% B) # Equation 4.20, vrt. fsr McDonald
  
  fsr <- factor.score.reliability(Lambda = B, Phi = diag(ncol(B)), Estimators = c("Regression", "Bartlett", "McDonald"))
  Regression <- fsr$Regression
  Bartlett <- fsr$Bartlett
  McDonald <- fsr$McDonald
}  
  if (sum) {
    WW <- matrix(rep(1,nrow(B)),nrow=nrow(B),ncol=1)
    a <- t(WW) %*% B %*% O %*% t(B) %*% WW
    b <- t(WW) %*% R %*% WW
    c <- diag(t(WW) %*% diag(diag(R-(B %*% t(B)))) %*% WW)
    reli <- rbind(reli,c(1/(1+c/a),a/b))
    rownames(reli) <- c(colnames(B),"Sum")
    oa <- c(oa,diag(t(WW) %*% RU %*% WW))
    ob <- c(ob,diag(t(WW) %*% R %*% WW))
    ab <- c(ab,sum(diag(S)))
    aa <- c(aa,t(WW) %*% S %*% WW)  

    if (scores) {    
    rhoScore <- c(rhoScore,NA)
    McDonald <- c(McDonald,NA)
    Regression <- c(Regression,NA)
    Bartlett <- c(Bartlett,NA)
    }
  }
  Alpha <- (p/(p-1))*(1-ab/aa)
  Omega <- c(oa/ob)
  
  reli <- cbind(reli,Alpha)
  if (!scores) { reli <- cbind(reli,Omega) }
  if (scores)  { 
    reli <- cbind(reli,rhoScore,McDonald,Regression,Bartlett)
  }

  return(reli)
}

source("factor_score_reliability.R")

colStds <- function(sample) { apply(sample,2,sd)}

rlbsimul <- function(M,SD,R,B,n,k,f) {
  sr1 <- NULL
  sr2 <- NULL
  for (i in 1:k) {
    S <- R * SD %*% t(SD)
    sample <- MASS::mvrnorm(n, mu = M, Sigma = S)
    sSD <- colStds(sample)
    sR <- cor(sample)
    A <- factanal(covmat=sR,factors=f)$loadings
    class(A) <- NULL
    UDV <- svd(t(A) %*% B)
    L <- UDV$u %*% UDV$v
    F <- A %*% L
    colnames(F) <- paste0("Image",1:ncol(F))
    r1 <- as.data.frame(reli(sR,F,M=sSD,sum=TRUE,score=FALSE))
    sW <- solve(sR) %*% F
    colnames(sW) <- paste0("Score",1:ncol(sW))
    r2 <- as.data.frame(reli(sR,F,sW,M=sSD))
    r1$lab <- rownames(r1)
    r1$sim <- i
    r2$lab <- rownames(r2)
    r2$sim <- i
    sr1 <- rbind(sr1,r1)
    sr2 <- rbind(sr2,r2)
  }
  return(list(images=sr1,scores=sr2))
}


simuexp <- function(expseed,simun,ssize_all,prilo_all) {
  
  ts <- NULL
  
  count <- 1
  pb = txtProgressBar(min = 1, max = length(ssize_all)*length(prilo_all), initial = 1) 
  
  
  for (i in 1:length(prilo_all)) {
    prilo <- prilo_all[i]  
    B <- genB(prilo)
    colnames(B) <- paste0("Image",1:ncol(B))
    rownames(B) <- paste0("Item",1:nrow(B))
    M <- rep(0,each=nrow(B))
    SD <- rep(1,each=nrow(B))
    R <- B %*% t(B) + (diag(nrow(B)) - diag(diag(B %*% t(B))))
    r1 <- as.data.frame(reli(R,B,M=SD,sum=TRUE,scores=FALSE))
    W <- solve(R) %*% B     # Weights for factor scores (orthogonal factors)
    colnames(W) <- paste0("Score",1:ncol(W))
    r2 <- as.data.frame(reli(R,B,W,M=SD))
    r1$lab <- rownames(r1)
    r2$lab <- rownames(r2)
    
    for (j in 1:length(ssize_all)) {
      ssize <- ssize_all[j]
      
      setTxtProgressBar(pb,count)
      ts$true[[count]] <- list(images=r1,scores=r2)
      set.seed(expseed+count)
      ts$simu[[count]] <- rlbsimul(M,SD,R,B,ssize,simun,ncol(B))
      count <- count +1
      
    }
  }
  
  return(ts)
  
}








