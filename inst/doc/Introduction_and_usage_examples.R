## ---- fig.width = 3.3, fig.height = 2.7----------------------------------
library(mdpeer)

n <- 100
p1 <- 10
p2 <- 90
p <- p1 + p2

# Define graph Adjacency matrix
A <- matrix(rep(0, p*p), nrow = p, ncol = p)
A[1:p1, 1:p1] <- 1
A[(p1+1):p, (p1+1):p] <- 1

# Compute graph Laplacian matrix 
L <- Adj2Lap(A)
L.norm <- L2L.normalized(L)

# Vizualize matrices
vizu.mat(A, title = "Adjacency matrix")
vizu.mat(L, title = "Laplacian matrix"); vizu.mat(L.norm, title = "Laplacian matrix (normalized)")

## ------------------------------------------------------------------------
set.seed(1234)
n <- 200 
p1 <- 10
p2 <- 90
p <- p1 + p2
A <- matrix(rep(0, p*p), nrow = p, ncol = p)
A[1:p1, 1:p1] <- 1
A[(p1+1):p, (p1+1):p] <- 1
L <- Adj2Lap(A)
L.norm <- L2L.normalized(L)
Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
b.true<- c(rep(1, p1), rep(0, p2))
beta.true <- runif(3)
intercept <- 0
eta <- intercept + Z %*% b.true 
R2 <- 0.5 # assumed variance explained 
sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
error <- rnorm(n, sd = sd.eps)
Y <- eta + error

## ---- message=FALSE------------------------------------------------------
?RidgePEER

RidgePEER.fit <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = NULL)

# Optimal lambda regularization parameter values
c(RidgePEER.fit$lambda.Q, RidgePEER.fit$lambda.R)

## ---- fig.width = 3.3, fig.height = 3.3----------------------------------
# Intercept estimate (the only non-penalized coefficient in this setting)
RidgePEER.fit$beta.est

# Compare true b estimates and RidgePEER estimates 
b.est.RidgePEER <- RidgePEER.fit$b.est
plot.y.lim <- range(c(b.true, b.est.RidgePEER))
par(cex = 0.7)
plot(b.true, main = "b.true", ylab = "", xlab = "", ylim = plot.y.lim); plot(b.est.RidgePEER, main = "b.est.RidgePEER", ylab = "", xlab = "", col = "blue", ylim = plot.y.lim)

# b estimation MSE 
mean((b.true - b.est.RidgePEER)^2)

## ------------------------------------------------------------------------
set.seed(1234)
n <- 200 
p1 <- 10
p2 <- 90
p <- p1 + p2
A <- matrix(rep(0, p*p), nrow = p, ncol = p)
A[1:p1, 1:p1] <- 1
A[(p1+1):p, (p1+1):p] <- 1
L <- Adj2Lap(A)
L.norm <- L2L.normalized(L)
Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
b.true<- c(rep(1, p1), rep(0, p2))
X <- matrix(rnorm(n*3), nrow = n, ncol = 3)
beta.true <- runif(3)
intercept <- 0
eta <- intercept + Z %*% b.true + X %*% beta.true
R2 <- 0.5 # assumed variance explained 
sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
error <- rnorm(n, sd = sd.eps)
Y <- eta + error

## ---- message=FALSE------------------------------------------------------
RidgePEER.fit <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X)

# Optimal lambda regularization parameter values
c(RidgePEER.fit$lambda.Q, RidgePEER.fit$lambda.R)

## ------------------------------------------------------------------------
# Intercept and 3 covariates estimates 
RidgePEER.fit$beta.est

## ---- fig.width = 3.3, fig.height = 3.3----------------------------------
# Compare true b estimates and RidgePEER estimates 
b.est.RidgePEER <- RidgePEER.fit$b.est
plot.y.lim <- range(c(b.true, b.est.RidgePEER))
par(cex = 0.7)
plot(b.true, main = "b.true", ylab = "", xlab = "", ylim = plot.y.lim); plot(b.est.RidgePEER, main = "b.est.RidgePEER", ylab = "", xlab = "", col = "blue", ylim = plot.y.lim)

# b estimation MSE 
mean((b.true - b.est.RidgePEER)^2)

## ------------------------------------------------------------------------
set.seed(1234)
n <- 200 
p1 <- 10
p2 <- 90
p <- p1 + p2
A <- matrix(rep(0, p*p), nrow = p, ncol = p)
A[1:p1, 1:p1] <- 1
A[(p1+1):p, (p1+1):p] <- 1
L <- Adj2Lap(A)
L.norm <- L2L.normalized(L)
Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
b.true <- c(rep(1, p1), rep(0, p2))
beta.true <- runif(3)
intercept <- 0
eta <- intercept + Z %*% b.true
R2 <- 0.15 # assumed variance explained 
sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
error <- rnorm(n, sd = sd.eps)
Y <- eta + error

## ---- fig.width = 3.3, fig.height = 2.7, echo = FALSE--------------------
par(cex = 0.7)
vizu.mat(A, title = "Adjacency matrix", base_size = 8); plot(b.true, main = "b.true", ylab = "", xlab = "")

## ---- message=FALSE------------------------------------------------------
RidgePEER.fit <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X)
PEER.fit      <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X, add.Ridge = FALSE)
Ridge.fit     <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X, add.PEER = FALSE)

## ---- echo=FALSE---------------------------------------------------------
par(cex = 0.7)
plot(RidgePEER.fit$b.est, main = "b coefficients \n(graph-orig. + Ridge-orig. penalty)", ylab = "", xlab = "")
plot(PEER.fit$b.est, main = "b coefficients \n(graph-originated penalty only)", ylab = "", xlab = "")
plot(Ridge.fit$b.est, main = "b coefficients \n(Ridge-originated penalty only)", ylab = "", xlab = "")

## ------------------------------------------------------------------------
# b coeffcient estimates MSE
RidgePEER.b.MSE <- mean((RidgePEER.fit$b.est - b.true)^2)
PEER.b.MSE      <- mean((PEER.fit$b.est - b.true)^2)
Ridge.b.MSE     <- mean((Ridge.fit$b.est - b.true)^2)

# MSE
MSE.vec <- c(RidgePEER.b.MSE, PEER.b.MSE, Ridge.b.MSE)
names(MSE.vec) <- c("RidgePEER", "PEER", "Ridge")
round(MSE.vec, 4)

# MSE: % of RidgePEER
round(MSE.vec*(1/MSE.vec[1]), 4)

## ------------------------------------------------------------------------
set.seed(1234)
n <- 200 
p1 <- 10
p <- p1*10
A <- matrix(rep(0, p*p), nrow = p, ncol = p)
A.pos <- as.logical(rep(c(rep(1, 10), rep(0, 10)), 5))
A[A.pos, A.pos] <- 1
L <- Adj2Lap(A)
L.norm <- L2L.normalized(L)
Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
b.true <- as.numeric(c(A.pos[6:100], A.pos[1:5]))
X <- matrix(rnorm(n*3), nrow = n, ncol = 3)
beta.true <- runif(3)
intercept <- 0
eta <- intercept + Z %*% b.true+ X %*% beta.true
R2 <- 0.15 # assumed variance explained 
sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
error <- rnorm(n, sd = sd.eps)
Y <- eta + error

## ---- fig.width = 3.3, fig.height = 2.7, echo = FALSE--------------------
par(cex = 0.7)
vizu.mat(A, title = "Adjacency matrix", base_size = 8); plot(b.true, main = "b.true", ylab = "", xlab = "")

## ---- message=FALSE------------------------------------------------------
RidgePEER.fit <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X)
PEER.fit      <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X, add.Ridge = FALSE)
Ridge.fit     <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X, add.PEER = FALSE)

## ---- echo=FALSE---------------------------------------------------------
par(cex = 0.7)
plot(RidgePEER.fit$b.est, main = "b coefficients \n(graph-orig. + Ridge-orig. penalty)", ylab = "", xlab = "")
plot(PEER.fit$b.est, main = "b coefficients \n(graph-originated penalty only)", ylab = "", xlab = "")
plot(Ridge.fit$b.est, main = "b coefficients \n(Ridge-originated penalty only)", ylab = "", xlab = "")

## ------------------------------------------------------------------------
# b coeffcient estimates MSE
RidgePEER.b.MSE <- mean((RidgePEER.fit$b.est - b.true)^2)
PEER.b.MSE      <- mean((PEER.fit$b.est - b.true)^2)
Ridge.b.MSE     <- mean((Ridge.fit$b.est - b.true)^2)

# MSE
MSE.vec <- c(RidgePEER.b.MSE, PEER.b.MSE, Ridge.b.MSE)
names(MSE.vec) <- c("RidgePEER", "PEER", "Ridge")
round(MSE.vec, 4)

# MSE: % of RidgePEER
round(MSE.vec*(1/MSE.vec[1]), 4)

