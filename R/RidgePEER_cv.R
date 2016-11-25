# rPEER function data objects consistency 
# 
# Validate rPEER function data objects dimensions consistency 
# 
# @param Q graph Laplacian matrix
# @param y response values matrix (n x 1)
# @param Z design matrix (n x p) modeled as random effects variables 
# @param X design matrix (n x k) modeled as fixed effects variables 
# @return NULL
# 
validate.RidgePEER.cv.data.objects <- function(Q, y, Z, X){
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(Q)[2] == dim(Z)[2])) stop("dim(Q)[2] != dim(Z)[2]")
  if (!is.null(X)) if(!(dim(y)[1] == dim(X)[1])) stop("dim(y)[1] != dim(X)[1]")
}



#' @title Graph-constrained estimation with regulariazation parameter selection via Cross-Validation
#' 
#' @description 
#' Perform penalized estimation where penalty term is a linear combination of 
#' graph-originated penalty and Ridge-originated penalty terms. Corresponding regularization
#' parameter values \eqn{\lambda} (\eqn{\lambda_Q, \lambda_R}) are selected via K-fold Cross-Validation
#' procedure. 
#' 
#' Assumes MSE of response variable values prediction to be minimized in the Cross-Validation procedure.
#' 
#' @param Q graph-originated penalty matrix \eqn{(p \times p)}, typically: a graph Laplacian matrix
#' @param y response values matrix \eqn{(n \times 1)}
#' @param Z design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling)
#' @param X design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling)
#' @param Q.normalized logical whether or not to use normalized version of a graph-originated penalty matrix
#' @param r.const scalar value of regularization parameter for a Ridge term by which matrix \code{Q} is corrected 
#' (used in situations when current value of \code{lambda.2} Cross-Validation parameter grid is less than it)
#' @param lambda.Q.grid vector of \code{lambda.Q} values to iterate over in Cross-Validation procedure 
#' @param lambda.2.grid vector of \code{lambda.2} values to iterate over in Cross-Validation procedure 
#' @param K number of Cross-Validation procedure folds
#' @param verbose logical whether or not message out information from function execution
#' 
#' @return 
#' \item{b.est}{vector with estimated values of \eqn{b} coefficients}
#' \item{beta.est}{vector with estimated values of \eqn{\beta} coefficients}
#' \item{ERR.mat.cv}{matrix of MSE errors obtained in Cross-Validation for each parameter grids combination considered}
#' \item{lambda.Q.cv}{value of \eqn{\lambda_Q} regularization parameter selected via Cross-Validation}
#' \item{lambda.R.cv}{\code{lambda.Q.cv} * \code{lambda.2.cv} value}
#' \item{lambda.2.cv}{value of \eqn{\lambda_2} regularization parameter selected via Cross-Validation}
#' 
#' @examples
#' # Example 1. 
#' 
#' set.seed(1234)
#' # graph Adjacency matrix
#' n <- 200 
#' p1 <- 10
#' p2 <- 90
#' p <- p1 + p2
#' A <- matrix(rep(0, p*p), nrow = p, ncol = p)
#' A[1:p1, 1:p1] <- 1
#' A[(p1+1):p, (p1+1):p] <- 1
#' # graph Laplacian matrix
#' L <- Adj2Lap(A)
#' L.norm <- L2L.normalized(L)
#' # Z design matrix 
#' Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' # True b coefficients 
#' b.true<- c(rep(1, p1), rep(0, p2))
#' beta.true <- runif(3)
#' intercept <- 0
#' eta <- intercept + Z %*% b.true 
#' R2 <- 0.5 # assumed variance explained 
#' sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
#' error <- rnorm(n, sd = sd.eps)
#' # Y observed 
#' Y <- eta + error
#' 
#' # model estimation
#' RidgePEER.cv.fit <- RidgePEER.cv(Q = L.norm, y = Y, Z = Z, X = NULL,
#'     lambda.Q.grid = exp(seq(-6, 6, length.out = 2)),
#'     lambda.2.grid = exp(seq(-6, 6, length.out = 2)))
#' 
#' # b coefficient estimates
#' RidgePEER.cv.fit$b.est
#' # regularization parameters values
#' RidgePEER.cv.fit$lambda.Q
#' RidgePEER.cv.fit$lambda.2 
#' RidgePEER.cv.fit$lambda.R # lambda.R = lambda.Q * lambda.2
#' 
#' 
#' # Example 2.: model with non-penalized covariates 
#' 
#' # X design matrix (covariates which are not to be penalized)
#' X <- matrix(rnorm(n*3), nrow = n, ncol = 3)
#' beta.true <- runif(3)
#' intercept <- 0
#' eta <- intercept + Z %*% b.true + X %*% beta.true
#' R2 <- 0.5 # assumed variance explained 
#' sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
#' error <- rnorm(n, sd = sd.eps)
#' # Y observed 
#' Y <- eta + error
#' 
#' # model estimation
#' RidgePEER.cv.fit <- RidgePEER.cv(Q = L.norm, y = Y, Z = Z, X = X,
#'     lambda.Q.grid = exp(seq(-6, 6, length.out = 2)),
#'     lambda.2.grid = exp(seq(-6, 6, length.out = 2)))
#' 
#' # b coefficient estimates
#' RidgePEER.cv.fit$b.est
#' # beta coefficient estimates
#' RidgePEER.cv.fit$beta.est
#' # regularization parameters values
#' RidgePEER.cv.fit$lambda.Q
#' RidgePEER.cv.fit$lambda.2 
#' RidgePEER.cv.fit$lambda.R # lambda.R = lambda.Q * lambda.2
#' 
#' 
#' @import glmnet
#' @export
#' 
#'
RidgePEER.cv <- function(Q, y, Z, X = NULL,
                         Q.normalized = FALSE,
                         r.const = 0.001, 
                         lambda.Q.grid = exp(seq(-6, 6, length.out = 10)),
                         lambda.2.grid = exp(seq(-6, 6, length.out = 10)),
                         K = 10, 
                         verbose = FALSE){
  
  # Check for data objects dimensions consistency 
  validate.RidgePEER.cv.data.objects(Q, y, Z, X)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  pcov <- ifelse(is.null(X), 0, dim(X)[2])
  
  # Transform data objects 
  if (Q.normalized) Q <- L2L.normalized(Q)
  
  # Regress out fixed effect coefficients (keep suffix ".wave" even if there is no fixed effect coefficients)
  if (pcov > 0){
    X.hash <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
    y.wave <- X.hash %*% y
    Z.wave <- X.hash %*% Z
  } else {
    y.wave <- y
    Z.wave <- Z
  }
  
  # Cross-validation
  ERR.mat <- matrix(NA, nrow = length(lambda.Q.grid), ncol = length(lambda.2.grid))
  
  for(iQ in 1:length(lambda.Q.grid)){
    for(i2 in 1:length(lambda.2.grid)){
      # iQ <- 5 ; i2 <- 5
      
      # Regularization parameters of current interation
      lQ.tmp <- lambda.Q.grid[iQ]
      l2.tmp <- lambda.2.grid[i2]
      if (l2.tmp < r.const){
        l2.tmp <- r.const
      }
      if (verbose){
        print("RidgePEER.cv:")
        print(paste0("   lQ: ", lQ.tmp))
        print(paste0("   l2: ", l2.tmp))
      }
      
      tryCatch({
        # Define LMM random effects matrix lme.Z
        Q_tilde <-  lQ.tmp * Q + l2.tmp * diag(p)
        Q_tilde.chol <- chol(Q_tilde)
        Q_tilde.chol.inv <- solve(Q_tilde.chol)
        lme.Z.wave <- Z.wave %*% Q_tilde.chol.inv
        
        # Cross-validation for a Ridge problem (with lambda set to 1 to reflect the previous situation)
        # We divide by n to reflect the penalty used in glment Ridge, see:
        # http://stackoverflow.com/questions/39863367/ridge-regression-with-glmnet-gives-different-coefficients-than-what-i-compute
        cv.lambda <- c(1, 1.0001) * (1/n)
        cv.fit <- cv.glmnet(x = lme.Z.wave, y = y.wave, alpha = 0, lambda = cv.lambda, nfolds = K, standardize = FALSE, intercept = FALSE)
        mse.cv.lambda <- cv.fit$cvm[cv.fit$lambda == (1/n)]
        ERR.mat[iQ, i2] <- mse.cv.lambda
        
      }, error = function(e) {
        if (verbose){
          message("RidgePEER.cv: error for:")
          message(paste0("lQ: ", lQ.tmp))
          message(paste0("l2: ", l2.tmp))
        }
      })
    }
  }
  
  # If all cross-validations were unsuccessful
  if (all(is.na(ERR.mat))){
    message("RidgePEER.cv: error - all cross-validations were unsuccessful. Returning NULLs")
    res.list <- list(b.est = NULL, beta.est = NULL, ERR.mat.cv = NULL, lambda.Q.cv = NULL, lambda.2.cv = NULL)
    return(res.list)
  }
  
  # Best cross-validation regularization parameters
  cv.best.idx <- which(ERR.mat == min(ERR.mat, na.rm = TRUE), arr.ind = TRUE)
  lambda.Q.cv <- lambda.Q.grid[cv.best.idx[1]]
  lambda.2.cv <- lambda.2.grid[cv.best.idx[2]]
  if (lambda.2.cv < r.const){
    lambda.2.cv <- r.const
  }
  
  # Compute coefficient estimates
  Q.cv <-  lambda.Q.cv * Q + lambda.2.cv * diag(p)
  b.est <- as.vector(solve(t(Z.wave) %*% Z.wave + Q.cv) %*% t(Z.wave) %*% y.wave)
  if (pcov > 0){
    beta.est <- as.vector(solve(t(X) %*% X) %*% t(X) %*% (y - Z %*% b.est))
  } else {
    beta.est <- c() 
  }
  
  res.list <- list(b.est = b.est, 
                   beta.est = beta.est, 
                   ERR.mat.cv = ERR.mat, 
                   lambda.Q.cv = lambda.Q.cv, 
                   lambda.R.cv = lambda.Q.cv * lambda.2.cv,
                   lambda.2.cv = lambda.2.cv)
  return(res.list)
}