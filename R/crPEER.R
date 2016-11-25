# Validate crPEER function data objects dimensions consistency 
# 
# @param Q graph Laplacian matrix
# @param y response values matrix (n x 1)
# @param Z design matrix (n x p) modeled as random effects variables 
# @param X design matrix (n x k) modeled as fixed effects variables 
# @return NULL
# 
validate.crPEER.data.objects <- function(Q, y, Z, X){
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(Q)[2] == dim(Z)[2])) stop("dim(Q)[2] != dim(Z)[2]")
  if (!is.null(X)) if(!(dim(y)[1] == dim(X)[1])) stop("dim(y)[1] != dim(X)[1]")
}



#' Graph-constrained regression with additional Ridge term to handle the non-invertibility of 
#' a graph-originated penalty matrix
#' 
#' @param Q graph-originated penalty matrix \eqn{(p \times p)}, typically: a graph Laplacian matrix
#' @param y response values matrix \eqn{(n \times 1)}
#' @param Z design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling)
#' @param X design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling)
#' @param lambda.2 scalar value of regularization parameter for additional Ridge term by which matrix \code{Q} is corrected 
#' @param Q.normalized logical whether or not to use normalized version of a graph-originated penalty matrix
#' @param compute.Grace logical whether or not to compute \href{https://arxiv.org/abs/1506.08339}{Grace test} 
#' (a significance test for graph-constrained estimation) results 
#' 
#' @return 
#' \item{b.est}{vector with estimated values of \eqn{b} coefficients}
#' \item{beta.est}{vector with estimated values of \eqn{\beta} coefficients}
#' \item{lambda.Q}{value of \eqn{\lambda_Q} regularization parameter}
#' \item{lambda.R}{\code{lambda.Q} * \code{lambda.2} value}
#' \item{lambda.2}{assumed "Ridge correction" \code{lambda.2} fixed regularization parameter value}
#' \item{grace.test.res}{output from the Grace significance test for graph-constrained estimation computation}
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
#' crPEER.fit <- crPEER(Q = L.norm, y = Y, Z = Z, X = NULL)
#' 
#' # b coefficient estimates
#' crPEER.fit$b.est
#' # regularization parameters values
#' crPEER.fit$lambda.Q
#' crPEER.fit$lambda.2 # assumed to be a fixed constant 
#' crPEER.fit$lambda.R # lambda.R = lambda.Q * lambda.2
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
#' crPEER.fit <- crPEER(Q = L.norm, y = Y, Z = Z, X = X)
#' 
#' # b coefficient estimates
#' crPEER.fit$b.est
#' # beta coefficient estimates
#' crPEER.fit$beta.est
#' # regularization parameters values
#' crPEER.fit$lambda.Q
#' crPEER.fit$lambda.2 # assumed to be a fixed constant 
#' crPEER.fit$lambda.R # lambda.R = lambda.Q * lambda.2
#' 
#' @references 
#' Brumback, B. A., Ruppert, D., Wand, M. P., Comment on 'Variable selection and function estimation in
#' additive nonparametric regression using a data-based prior'. Journal of the American Statistical
#' Association (1999): 94, 794–797.
#' 
#' Karas, M., Brzyski, D., Randolph, T., Harezlak, D. Brain connectivity-informed regularization methods for regression. Paper in progress, to be submited as an invited paper on CCNS for a special issue of  Statistics in Biosciences by Nov 30, 2016 (reference will be updated).
#'
#' Li, C., Li, H., Network-constrained regularization and variable selection for analysis of genomic data. 
#' Bioinformatics (2008): 24(9), 1175-1182.
#' 
#' Li, C., Li, H., Variable selection and regression analysis for graph-structured covariates with an application to genomics. 
#' The Annals of Applied Statistics (2010): 4(3), 1498–1516.
#' 
#' Randolph, T., Harezlak, J., Feng, Z., Structured penalties for functional linear models—partially empirical 
#' eigenvectors for regression. The Electronic Journal of Statistics (2012): 6, 323-353.
#' 
#' @import nlme
#' @export
#' 
crPEER <- function(Q, y, Z, X = NULL, lambda.2 = 0.001, Q.normalized = FALSE, 
                   compute.Grace = FALSE){
  
  # Check for data objects dimensions consistency 
  validate.crPEER.data.objects(Q, y, Z, X)
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
  
  # Define LMM random effects matrix lme.Z
  Q_tilde <-  Q + lambda.2 * diag(p)
  Q_tilde.chol <- chol(Q_tilde)
  Q_tilde.chol.inv <- solve(Q_tilde.chol)
  lme.Z.wave <- Z.wave %*% Q_tilde.chol.inv
  
  # Fit LMM
  id.bd1 <- factor(rep(1, n))
  lme.fit <- lme(fixed = y.wave ~ 1, random = list(id.bd1 = pdIdent(~ lme.Z.wave - 1)), method = "REML") 
  
  # LMM-originated lambda parameter 
  sigma.eps <- lme.fit$sigma 
  sigma.u <- as.numeric(as.matrix(VarCorr(lme.fit))[1, 2])
  lambda.Q <- (sigma.eps^2)/(sigma.u^2)
  # Define lambda.R 
  lambda.R <- lambda.Q * lambda.2
  
  # Compute Grace Test results
  if (!compute.Grace){
    grace.test.res <- NULL
  } else {
    grace.test.res <- tryCatch({
      Q.opt <- lambda.Q * Q + lambda.R * diag(p)
      gt.res <- grace.test.C(Y = y.wave, X = Z.wave, L = Q.opt, lambda.L = 1, lambda.2 = 0, normalize.L = F,
                             normalize.X  = FALSE, center.Y = FALSE)
      list(b.est = gt.res$beta, b.est.pval = as.vector(gt.res$pvalue), Q.opt = Q.opt)
    }, error = function(e) {
      msg <- "ERROR OCURRED WHEN PROCEEDING WITH grace.test"
      message(msg); message(e)
      msg
    })  
  }
  
  # Compute coefficient estimates
  b.est <- as.vector(solve(t(Z.wave) %*% Z.wave + lambda.Q * Q_tilde) %*% t(Z.wave) %*% y.wave)
  if (pcov > 0){
    beta.est <- as.vector(solve(t(X) %*% X) %*% t(X) %*% (y - Z %*% b.est))
  } else {
    beta.est <- c() 
  }
  
  res.list <- list(b.est = b.est, 
                   beta.est = beta.est, 
                   lambda.Q = lambda.Q,
                   lambda.2 = lambda.2,
                   lambda.R = lambda.R,
                   grace.test.res = grace.test.res)
  return(res.list)
}