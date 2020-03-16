###############################################################################
#
# Original code from: André Beauducel, Christopher Harms & Norbert Hilger (2016).
# Reliability Estimates for Three Factor Score Estimators. International Journal
# of Statistics and Probability 5 (6), 94-107. doi: 10.5539/ijsp.v5n6p94
#
# They present an implementation of equation 4 in: Norman Cliff (1988). The
# Eigenvalues-Greater-Than-One Rule and the Reliability of Components.
# Psychological Bulletin 103(2), 276-279. The basic idea is that only the
# parameters of the factor model are necessary in order to calculate the
# reliabilities, when the hypothetical item set is equivalent.
#
# R code modified to work with one-dimensional case by
#   Reijo Sund (ORCID: 0000-0002-6268-8117), https://connect.uef.fi/en/reijo.sund
#
###############################################################################

##' This function computes and returns reliability estimates for three commonly used
##' Factor Score Estimators in Factor Analyses.
##'
##' Explanations of the algebraic formulas are presented in the paper. An
##' example how to use the function is given at the end of the file.
##'
##' @title Function for calculating reliability estimates for factor score Estimators
##' @param Lambda a \code{matrix} containing the loadings of items on the factors
##' @param Phi a \code{matrix} containing the factor intercorrelations
##' @param Estimators a \code{vector} to select the Estimators for which the reliability
##' 	estimates should be calculated. Available values: \code{Regression}, \code{Bartlett},
##' 	\code{McDonald}
##' @return Returns a two-dimensional list containing the reliability estimates for each
##' 	factor. Depending on the \code{Estimators} parameter, the list contains the values
##'		only for the selected Estimators.
##' @export
##' @author André Beauducel (\email{beauducel@uni-bonn.de})
##' @author Christopher Harms (\email{christopher.harms@uni-bonn.de})
##' @author Norbert Hilger (\email{nhilger@uni-bonn.de})
##'

factor.score.reliability <- function(Lambda, Phi, Estimators=c("Regression", "Bartlett", "McDonald")) {
  # Helper functions for frequently used matrix operations
  Mdiag <- function(x) return(diag(diag(x),ncol=ncol(x)))
  inv <- function(x) return(solve(x))
  
  # If a 'loadings' class is provided for lambda, we can easily convert it
  if (is(Lambda, "loadings"))
#    Lambda <- Lambda[,]
    class(Lambda) <- NULL
  
  # Perform several validity checks of the provided arguments
  if (any(missing(Lambda), missing(Phi), is.null(Lambda), is.null(Phi)))
    stop("Missing argument(s).")
  if (any(nrow(Phi) == 0, nrow(Lambda) == 0, ncol(Phi) == 0, ncol(Lambda) == 0))
    stop("Some diemension(s) of Phi or Lambda seem to be empty.")
  if (nrow(Phi) != ncol(Phi))
    stop("Phi has to be a q x q matrix.")
  if (ncol(Lambda) != nrow(Phi))
    stop("Phi and Lambda have a different count of factors.")
  if (any(round(min(Phi)) < 0, round(max(Phi)) > 1))
    stop("Phi contains invalid values (outside [0; 1]).")
  Estimators.Allowed <- c("Regression", "Bartlett", "McDonald")
  if (is.null(Estimators)) {
    message("No 'Estimators' defined, use 'Regression' as default.")
    Estimators <- c("Regression")
  }
  Estimators <- match.arg(Estimators, Estimators.Allowed, several.ok = TRUE)
  
  # Regenerate covariance matrix from factor loadings matrix
  Sigma <- (Lambda %*% Phi %*% t(Lambda))
  Sigma <- Sigma - Mdiag(Sigma) + diag(nrow(Lambda))
  
  # Calculate uniqueness/error of items
  Psi <- Mdiag(Sigma - Lambda %*% Phi %*% t(Lambda))^0.5
  
  ret <- list()
  if ("Regression" %in% Estimators) {
    # Reliability of Thurstone's Regression Factor Score Estimators
    # cf. Equation 6 in manuscript
    Rtt.Regression <- inv( Mdiag( Phi %*% t(Lambda) %*% inv(Sigma) %*% Lambda %*% Phi ) )^0.5 %*%
      Mdiag( Phi %*% t(Lambda) %*% inv(Sigma) %*% Lambda %*% Phi %*% 
               t(Lambda) %*% inv(Sigma) %*% Lambda %*% Phi) %*%
      inv( Mdiag( Phi %*% t(Lambda) %*% inv(Sigma) %*% Lambda %*% Phi ) )^0.5
    ret$Regression <- diag(Rtt.Regression)
  }
  if ("Bartlett" %in% Estimators) {
    # Reliability of Bartlett's Factor Score Estimators
    # cf. Equation 11 in manuscript
    Rtt.Bartlett <- inv( Mdiag( inv(t(Lambda) %*% inv(Psi)^2 %*% Lambda) + Phi ) )
    ret$Bartlett <- diag(Rtt.Bartlett)
  }
  if ("McDonald" %in% Estimators) {
    # Reliability of McDonald's correlation preserving factor score Estimators
    # cf. Equation 12 in manuscript
    Decomp <- svd(Phi)
    N <- Decomp$u %*% abs(diag(Decomp$d))^0.5
    sub.term <- t(N) %*% t(Lambda) %*% inv(Psi)^2 %*% Sigma %*% inv(Psi)^2 %*% Lambda %*% N
    Decomp <- svd(sub.term)
    d <- Decomp$d
    if (length(d)==1) { d <- as.matrix(Decomp$d) }
    sub.term <- Decomp$u %*% (diag(d)^0.5) %*% t(Decomp$u)
    
    Rtt.McDonald <- Mdiag( inv(sub.term) %*% t(N) %*% t(Lambda) %*% inv(Psi)^2 %*% Lambda %*%
                             Phi %*% t(Lambda) %*% inv(Psi)^2 %*% Lambda %*% N %*% inv(sub.term))
    ret$McDonald <- diag(Rtt.McDonald)
  }
  
  # Return reliabilities as list, so it can be accessed via e.g. factor.score.reliability(L, P)$Regression
  return(ret)
}

# # Example 1:
# Loadings <- matrix(c(
# 	 0.50,-0.10, 0.10,
# 	 0.50, 0.10, 0.10,
# 	 0.50, 0.10,-0.10,
# 	-0.10, 0.50, 0.15,
# 	 0.15, 0.50, 0.10,
# 	-0.15, 0.50, 0.10,
# 	 0.10, 0.10, 0.60,
# 	 0.10,-0.10, 0.60,
# 	 0.10, 0.10, 0.60
#  	),
# 	nrow=9, ncol=3,
# 	byrow=TRUE)
# InterCorr <- matrix(c(
# 	1.00, 0.30, 0.20,
#  	0.30, 1.00, 0.10,
#  	0.20, 0.10, 1.00
# 	),
# 	nrow=3, ncol=3,
# 	byrow=TRUE)
# 
# reliabilities <- factor.score.reliability(Lambda = Loadings, Phi = InterCorr, Estimators = c("Regression", "Bartlett", "McDonald"))
# lapply(reliabilities, round, 6)
# 
# factor.score.reliability(Lambda = B, Phi = diag(ncol(B)), Estimators = c("Regression", "Bartlett", "McDonald"))
# 
################################################## 
#
# Testing with reli, not confirmed to work / RS
# Bex <- Loadings
# Rex <- B %*% t(B) + (diag(nrow(B)) - diag(diag(B %*% t(B))))
# Oex <- InterCorr
# Wex <- solve(R) %*% B %*% O    # Weights for factor scores (correlated factors)
# 
# reli(Rex,Bex,Wex,Oex)
# 
# reli(Rex,Bex,O=Oex)





