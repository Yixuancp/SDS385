# SDS 385 Excercise 1 
# Generalized linear models

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


### (B) Steepest Gradient Method
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
liketrace <- likelihood(X, y, b)


# Iteration
for (step in 1:50000) {
  stepsize <- 0.02
  b <- b - stepsize * ll.grad(X, y, b) 
  lltrace <- c(lltrace,ll(X, y, b)) 
}
plot(lltrace)


### (D) Newton's method

# Define Hessian matrix 
ll.hes <-  function(X, b) {
  crossprod(diag(as.vector(sqrt(w(X,b) * (1 - w(X,b))))) %*% X)
}

# Initiation
b <- coef(lm(y ~ 0 + X))
#b <- matrix(0,11,1)
btrace <- b
lltrace <- ll(X, y, b)

# Iteration
for (step in 1:10) {
  b <- b - solve(ll.hes(X,b)) %*% ll.grad(X, y, b)
  lltrace <- c(lltrace,ll(X, y, b)) 
}
plot(lltrace)








