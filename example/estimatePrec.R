library(methods)
library(EmpericalSparsePrecision)
library(Matrix)
N <- 200
d <- 100
a <- 0.5
i <- c( 1:d, 1:(d-1), 2:d )
j <- c( 1:d, 2:d,     1:(d-1))
x <- c( c(1,rep(1 + a^2,d-2),1), rep(-a,d - 1), rep(-a, d - 1))
Q <- sparseMatrix( i = i, j = j, x = x)
L <- t(chol(Q))
Sigma <- as.matrix(solve(Q))


# simulating data
X <- matrix(0,nrow=N, ncol = d)
for(i in 1:N)
  X[i,1:d]  <- as.matrix(solve(t(L), rnorm(d)))

Sigma_emp <-cov(X)


Q_emp <- CovToPrecision(Sigma_emp, Q)
print(paste("reg || Q - Sigma_emp^-1|| = ", norm(solve(Sigma_emp)-Q), sep=""))
print(paste("reg || Q - Q_emp|| = ", norm(Q_emp-Q), sep=""))
