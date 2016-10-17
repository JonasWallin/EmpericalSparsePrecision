library(Matrix)
#'
#'@title CovToPrecision
#'@name CovToPrecision
#'
#'@description
#' verify Q, Sigma is symmetric, Perm is a vector of length D
#' Estimate a sparse Precision matrix (inverse covariance matrix)
#' using the Emperical covariance matrix and the sparsity structure of the
#' precision matrix.
#'
#'
#'@param: Sigma       - (d x d) matrix the covariance matrix
#'@param: Q           - (d x d) symbolic matrix represnting graph structure
#'@param: Perm        - (d x 1) permutation list
#'@param: calc.reo    - (bool) if Perm is NULL and calc.reo is TRUE caclulate a reodering
#'@return a sparse (d x d) matrix
#'
#' @examples
#' library(Matrix)
#' d <- 10
#' a <- 0.9
#' i <- c( 1:d, 1:(d-1), 2:d )
#' j <- c( 1:d, 2:d,     1:(d-1))
#' x <- c( c(1,rep(1 + a^2,d-2),1), rep(-a,d - 1), rep(-a, d - 1))
#' Q <- sparseMatrix( i = i, j = j, x = x)
#' Sigma <- as.matrix(solve(Q))
#' Qback <- CovToPrecision(Sigma, Q)
#'
#' @export
#'
CovToPrecision <- function(Sigma, Q, Perm = NULL, calc.reo = TRUE)
{
  if(is.null(Perm) & calc.reo == TRUE)
    Perm <- reorder_cpp(Q)$Perm
  else if(is.null(Perm) & calc.reo == FALSE)
    Perm <- 1:dim(Sigma)[1]

  L <- CovToPrecisionChol(Sigma[Perm, Perm], Q[Perm, Perm])
  P <- L %*% Matrix::t(L)
  P[Perm, Perm] <- P
  return(P)
}

#'@title CovToPrecisionChol
#'
#'@param: Sigma - (d x d) matrix the covariance matrix
#'@param: Q     - (d x d) symbolic matrix represnting graph structure
#'@param: Lcol  - (d x 1) Lcol[[i]] list of non-zero values for row i
#'
#' @export
CovToPrecisionChol <- function(Sigma, Q)
{
  if(isSymmetric(Sigma)==FALSE)
    stop("Sigma must be symmetric")
  if(Matrix::isSymmetric(Q)==FALSE)
    stop("Q must be symmetric")

  #TODO: check Q sparse and symmetric
  Lsymbol <- symbolicCholeskyFactor_cpp(Q)
  L <- CovToPrecisionChol_cpp(Lsymbol, Sigma)
  return(L)
}

#'
#'@title CovToPrecisionCholR
#'
#'@param: Sigma - (d x d) matrix the covariance matrix
#'@param: Q     - (d x d) symbolic matrix represnting graph structure
#'@param: Lcol  - (d x 1) Lcol[[i]] list of non-zero values for row i
#'
CovToPrecisionChol_R <- function(Sigma, Q)
{
  #TODO: check Q sparse and symmetric
  Lsymbol <- symbolicCholeskyFactor_cpp(Q)
  d <- dim(Sigma)[1]
  Ddiag <- vector(mode="numeric", length=d)
  Tmat  <- sparseMatrix( i = 1:d, j = 1:d, x = rep(1,d))
  for(j in 1:(d-1)){
    index <- setdiff(which(Lsymbol[,j] >0),j)
    if(length(index)> 0){
      Sigma_LL <- Sigma[index, index]
      Sigma_jL <- Sigma[j, index]
      Tmat[index, j] <- - solve(Sigma_LL, Sigma_jL)
      Ddiag[j]       <- Sigma[j, j] + Sigma_jL%*% Tmat[index,j]
    }else{ Ddiag[j]  <- Sigma[j, j]}
    Ddiag[d]       <- Sigma[d, d]
  }
  L <- Tmat%*%sparseMatrix( i = 1:d, j = 1:d, x = 1/sqrt(Ddiag))
  return(L)
}
