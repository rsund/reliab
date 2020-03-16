###############################################################################
#
# R-replication of similar simulation studies as reported in 
#   Vehkalahti (2000). http://hdl.handle.net/10138/10570
#
# Author(s):
# Reijo Sund (ORCID: 0000-0002-6268-8117), https://connect.uef.fi/en/reijo.sund
#
###############################################################################

source("reli-fun.R")
simun <- 500
prilo_all <- seq(0.60,0.95,0.05)
ssize_all <- seq(40,300,20)
autosave <- TRUE

# Experimental setting I
genB <- function(prilo) {
  B <- matrix(rep(prilo,each=10))
  return(B)
  }
exp1 <- simuexp(expseed=1111,simun,ssize_all,prilo_all)
if (autosave) save(exp1,file="exp1.RData")

# Experimental setting II
genB <- function(prilo) {
  B <- matrix(rep(prilo,each=25))
  return(B)
}
exp2 <- simuexp(expseed=2222,simun,ssize_all,prilo_all)
if (autosave) save(exp2,file="exp2.RData")

# Experimental setting III
genB <- function(prilo) {
  B <- rbind(
    matrix(rep(prilo,each=10)),
    matrix(rep(0.3,each=5)),
    matrix(rep(-0.3,each=5))
  )
  return(B)
}
exp3 <- simuexp(expseed=3333,simun,ssize_all,prilo_all)
if (autosave) save(exp3,file="exp3.RData")

# Experimental setting IV
genB <- function(prilo) {
  B <- rbind(
    matrix(rep(cbind(prilo,-0.1),each=20),ncol=2),
    matrix(rep(cbind(0.2,prilo),each=10),ncol=2)
  )
  return(B)
}
start_time <- Sys.time()
exp4 <- simuexp(expseed=4444,simun,ssize_all,prilo_all)
end_time <- Sys.time()
end_time - start_time
if (autosave) save(exp4,file="exp4.RData")

# Experimental setting V
genB <- function(prilo) {
  B <- rbind(
    matrix(rep(cbind(prilo,0,0),each=20),ncol=3),
    matrix(rep(cbind(0,prilo,0),each=10),ncol=3),
    matrix(rep(cbind(0,0,prilo),each=5),ncol=3)
  )
  return(B)
}
exp5 <- simuexp(expseed=5555,simun,ssize_all,prilo_all)
if (autosave) save(exp5,file="exp5.RData")


