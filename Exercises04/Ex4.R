# SDS 385 Excercise 4

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

N <- dim(X)[1]
p <- dim(X)[2]


######## Improving SGD for logistic regression ########
w <- function(x, b) {
  1 / (1 + exp(- x %*% b))
}

ll <-  function(x, y, b) {
  -log(prod(w(x, b)^y) * prod((1 - w(x,b))^(1 - y)))
}


# Define the gradient function
ll.grad <-  function(x, y, b) {
  -x %*% (y - w(x,b))
}




#### (B) ####
  
# Initiation
b <- matrix(0,p,1)
lltrace <- ll(X, y, b)
G <- diag(0, p, p)


# Iteration
for (step in 1:10000) {
  s <- sample(1:N, 1)
  G <-  sqrt(G^2 + diag(as.numeric(ll.grad(X[s, ], y[s], b)^2)))
  b <- b - solve(G) %*% ll.grad(X[s, ], y[s], b) 
  lltrace <- c(lltrace,ll(X, y, b)) 
}
plot(lltrace, cex = 0.5)



