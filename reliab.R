###############################################################################
#
# R-replication of reliability analyses presented in Kimmo Vehkalahti (2000):
# RELIABILITY OF MEASUREMENT SCALES - Tarkkonen’s general method supersedes Cronbach’s alpha.
# Statistical Research Raports 17, The Finnish Statistical Society, Helsinki, 2000.
# http://hdl.handle.net/10138/10570
#
# Author(s):
# Reijo Sund (ORCID: 0000-0002-6268-8117), https://connect.uef.fi/en/reijo.sund
#
###############################################################################
source("reli-fun.R")

source("data.R")

# DECA: Measuring the physical capacity by decathlon scores
# Means, standard deviations, correlations (M, R)
# Rotated factor matrix (5 factors) and factor score coefficient matrix (B,W)
# Schemes 5.2 and 5.4 (Survo solutions with maximum likelihood estimation and VARIMAX-rotation)
print(decaR)
print(decaB)
print(decaM)
print(decaW)

# Reliabilities of the factor images and the unweighted sum of DECA
# Scheme 5.5
reli(decaR,decaB,M=decaM[,"stddev"],sum=TRUE,scores=FALSE)

# Reliabilities of the factor scores of DECA.
# Scheme 5.7
reli(decaR,decaB,decaW,M=decaM[,"stddev"])

# Example of Religiosity and Fatalism by Heise & Bohrnstedt (1970)
# Schemes 5.26 - 5.29
reli(CRF,FRF,M=MRF[,"stddev"],sum=TRUE,scores=FALSE)




# Defining the structure of the measurement model
# Scheme 6.1
B <- matrix(
  c(
    0.8,0.8,0.5,0.5,0.1,0.1,0.1,
    0.0,0.0,0.2,0.2,0.7,0.6,0.6
    ),
  nrow=7,
  ncol=2
  )
colnames(B) <- paste0("Image",1:ncol(B))
rownames(B) <- paste0("Item",1:nrow(B))

# Forming the matrices of means, standard deviations and correlations
# Scheme 6.2

M <- rep(0,nrow(B))
SD <- rep(1,nrow(B))
R <- B %*% t(B) + (diag(nrow(B)) - diag(diag(B %*% t(B))))


# Computing the true values of reliabilities of the factor images and the unweighted sum
# Scheme 6.4
reli(R,B,M=SD,sum=TRUE,scores=FALSE)

# Computing and displaying the factor score coefficients
# Scheme 6.5
W <- solve(R) %*% B     # Weights for factor scores (orthogonal factors)
colnames(W) <- paste0("Score",1:ncol(W))
print(W)

# Computing the true values of reliabilities of the factor scores
# Scheme 6.6
reli(R,B,W,M=SD) 

#W <- solve(R) %*% B %*% O    # Weights for factor scores (correlated factors)
#reli(R,B,W,O)


# Generating a multinormal random sample and computing the sample statistics
# Scheme 6.7
set.seed(1111) # Set the seed of the random number generator
n=100
S <- R * SD %*% t(SD)
sample <- MASS::mvrnorm(n, mu = M, Sigma = S)

sM <- colMeans(sample)
sSD <- colStds(sample)
sR <- cor(sample)

print(sM)
print(sSD)
print(sR)

# Computing the maximum likelihood factor analysis
# Scheme 6.8
A <- factanal(sample,factors=2)$loadings
class(A) <- NULL
print(A)

# Computing the symmetric transformation analysis
# Scheme 6.9
UDV <- svd(t(A) %*% B)
L <- UDV$u %*% UDV$v
res <- A %*% L - B
print(L)
print(res)

# Multiplying the factor matrix by the transformation matrix
# Scheme 6.10
F <- A %*% L
colnames(F) <- paste0("Image",1:ncol(F))
print(F)

#RF <- solve(t(L) %*% L)
#print(RF)

# Computing the reliabilities of the factor images and the unweighted sum from the sample
# Scheme 6.11
reli(sR,F,M=sSD,sum=TRUE,score=FALSE)

# Computing and displaying the factor score coefficients from the sample
# Scheme 6.12
sW <- solve(sR) %*% F     # Weights for factor scores (orthogonal factors)
colnames(sW) <- paste0("Score",1:ncol(sW))
print(sW)

# Computing the reliabilities of the factor scores from the sample
# Scheme 6.13
reli(sR,F,sW,M=sSD) 




# Making the central computations in a condensed form
# Scheme 6.14
# Defining the structure of the measurement model
B <- matrix(
  c(
    0.8,0.8,0.5,0.5,0.1,0.1,0.1,
    0.0,0.0,0.2,0.2,0.7,0.6,0.6
  ),
  nrow=7,
  ncol=2
)

source("prep-true.R")

# True values for factor images
print(t(r1))

# True values for factor scores
print(t(r2))


# Conducting a repeated simulation experiment
# Scheme 6.15
set.seed(1111) # Set the seed of the random number generator
n=100
sim <- rlbsimul(M,SD,R,B,n,50,2)


# True values of reliabilities and means of 50 sample estimates
# Table 6.3

# Factor Images
library(dplyr)

#True
print(t(r1))

# Simulated
sim$images %>%
  group_by(lab) %>%
  summarise_all(list(mean=mean,sd=sd)) %>%
  t()


#Factor Scores

#True
print(t(r2)) 

# Simulated
sim$scores %>%
  group_by(lab) %>%
  summarise_all(list(mean=mean,sd=sd)) %>%
  t()








