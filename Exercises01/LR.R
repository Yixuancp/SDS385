# SDS 385 Excercise 1 

library(microbenchmark)
library(Matrix)

### Linear regression ###
N <- 1000
p <- 50

X <- matrix(rnorm(N * p), nrow = N)
y <- matrix(rnorm(N), nrow = N)
W <- diag(1, N, N)

# Inversion Method
InvM <- function(X, y, W){ 
  beta.inv <- solve(crossprod(sqrt(W) %*% X)) %*% t(X) %*% W %*% y
  return(beta.inv)
}

# Cholesky Decomposition
CholD <- function(X, y, W){ 
  U <- chol(crossprod(sqrt(W) %*% X))
  L <- t(U)
  temp <- solve(L, t(X) %*% W %*% y)
  beta.cho <- solve(U, temp)
  return(beta.cho)
}

# Compare two methods
microbenchmark(
  InvM(X, y, W),
  CholD(X, y, W),
  times = 10
)

### Sparse Matrix ###
d <- 0.05
X.sparse <- matrix(rnorm(N * p), nrow = N)
mask <- matrix(rbinom(N * p, 1, d), nrow = N) 
X.sparse <- mask * X.sparse
W <- diag(1, N, N)

# Cholesky Decomposition for sparse matrix
CholD.sparse <- function(X, y, W){ 
  X <- Matrix(X, sparse = TRUE)
  U <- chol(crossprod(sqrt(W) %*% X))
  L <- t(U)
  temp <- solve(L, t(X) %*% W %*% y)
  beta.cho <- solve(U, temp)
  return(beta.cho)
}

microbenchmark(
  InvM(X.sparse, y, W),
  CholD(X.sparse, y, W),
  CholD.sparse(X.sparse, y, W),
  times = 10
)

