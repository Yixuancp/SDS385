# SDS 385 Excercise 2 
# Data from Ex1

rm(list=ls())
# Import and preprocess data
setwd("~/Box Sync/SDS385/Ex1")
wdbc <- read.csv("wdbc.csv", header=FALSE)

y = wdbc[,2]
X = as.matrix(wdbc[,-c(1,2)])
scrub = which(1:ncol(X) %% 3 == 0)
scrub = 11:30
X = X[,-scrub]
y<-ifelse(y=="B", 0, 1)
y<-as.numeric(y)
X <- scale(X)
glm1 = glm(y~X, family='binomial')
X <- cbind(1,X) # Add a column of 1s

N <- dim(X)[1]
p <- dim(X)[2]

w <- function(x, b) {
  1 / (1 + exp(- x %*% b))
}

ll <-  function(x, y, b) {
  -log(prod(w(x, b)^y) * prod((1 - w(x,b))^(1 - y)))
}

# Define the gradient function
ll.grad <-  function(x, y, b) {
  - x %*% (y - w(x,b))
}

# Constant stepsize
b <- matrix(0, p, 1)
lltrace <- ll(X, y, b)
stepsize <- 0.02
for (step in 1:50000) {
  s <- sample(1:N, 1)
  b <- b - stepsize * ll.grad(X[s, ], y[s], b) 
  lltrace <- c(lltrace, ll(X, y, b)) 
}
plot(lltrace, cex = 0.2)

# Robbins-Monro
bSGD <- matrix(0, p, 1) 
lltraceSGD <- ll(X, y, bSGD) # Track the loglikelihood
lltraceSGDPR <- ll(X, y, bSGD) # Track the loglikelihood using Polyak-Ruppert averaging
bMean <- bSGD # Average of estimates
alpha <- 0.9
C <- 2
for (step in 1:50000) {
  stepsize <- C * (step + 1)^(-alpha)
  s <- sample(1:N, 1)
  bSGD <- bSGD - stepsize * ll.grad(X[s, ], y[s], bSGD) 
  bMean <- (bMean * step + bSGD) / (step + 1)
  lltraceSGD <- c(lltraceSGD,ll(X, y, bSGD)) 
  lltraceSGDPR <- c(lltraceSGDPR,ll(X, y, bMean)) 
}
plot(lltraceSGD, cex = 0.2)
par(new = TRUE)
plot(lltraceSGDPR, cex = 0.2, col = "blue")


