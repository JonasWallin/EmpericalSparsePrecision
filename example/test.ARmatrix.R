library(methods)
library(EmpiricalSparsePrecision)
library(Matrix)
d <- 100
a <- 0.9
i <- c( 1:d, 1:(d-1), 2:d )
j <- c( 1:d, 2:d,     1:(d-1))
x <- c( c(1,rep(1 + a^2,d-2),1), rep(-a,d - 1), rep(-a, d - 1))
Q <- sparseMatrix( i = i, j = j, x = x)
Q[d-1, 2] <- .1
Q[2, d-1] <- .1
L <- t(chol(Q))
Sigma <- as.matrix(solve(Q))
Q_out <- CovToPrecision(Sigma, Q)
