# SDS 385 Excercise 3

rm(list=ls())
# Import and preprocess data
setwd("~/Box Sync/SDS385/Ex1")
wdbc <- read.csv("wdbc.csv", header=FALSE)
y = wdbc[ , 2]
X = as.matrix(wdbc[,-c(1,2)])
scrub = which(1:ncol(X) %% 3 == 0)
scrub = 11:30
X = X[,-scrub]
y <- ifelse(y == "B", 0, 1)
y <- as.numeric(y)
X <- scale(X)
glm1 = glm(y ~ X, family = 'binomial')
X <- cbind(1, X) # Add a column of 1s


######## Linear search ########

#### (A) ####




#### (B) Implementation ####
# Define the success probability
w <- function(x, b) {
  1 / (1 + exp(- x %*% b))
}

# Define the log likelihood function
ll <-  function(x, y, b) {
  -log(prod(w(x, b)^y) * prod((1 - w(x,b))^(1 - y)))
}

# Define the likelihood function
likelihood <-  function(x, y, b) {
  prod(w(x, b)^y) * prod((1 - w(x,b))^(1 - y))
}

# Define the gradient function
ll.grad <-  function(x, y, b) {
  -t(x) %*% (y - w(x,b))
}

# Intial value
b<- matrix(0, 11, 1)
btrace <- b
lltrace <- ll(X, y, b)
stepsizebar <- 0.1
rho <- 0.75

# Iteration
for (step in 1:50000) {
  stepsize <- stepsizebar
  for (s in 1:100){
  stepsize <- stepsize * rho
  bt <- b - stepsize * ll.grad(X, y, b) 
  if (ll(X, y, bt) < ll(X, y, b)) break
  }
  b <- bt
  lltrace <- c(lltrace,ll(X, y, b)) 
}
plot(lltrace, cex = 0.5)



######## Quasi-Newton's method (BFGS) ########

#### (A) ####



#### (B) Implementation ####
# Define Hessian matrix 
ll.hes <-  function(X, b) {
  crossprod(diag(as.vector(sqrt(w(X,b) * (1 - w(X,b))))) %*% X)
}

# Initiation
b <- matrix(0,11,1)
btrace <- b
lltrace <- ll(X, y, b)
B <- ll.hes(X,b)

# Iteration
for (step in 1:30) {
  sk <- - solve(B) %*% ll.grad(X, y, b)
  yk <- ll.grad(X, y, b + sk) - ll.grad(X, y, b)
  b <-  b + sk
  B <- B - crossprod(t(sk) %*% B) / as.numeric(t(sk) %*% B %*% sk) + tcrossprod(yk) / as.numeric(t(yk) %*% sk)
  lltrace <- c(lltrace,ll(X, y, b)) 
}
plot(lltrace, cex = 0.5)


