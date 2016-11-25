# vrPEER function data objects consistency 
# 
# Validate vrPEER function data objects dimensions consistency 
# 
# @param Q graph Laplacian matrix
# @param y response values matrix (n x 1)
# @param Z design matrix (n x p) modeled as random effects variables 
# @param X design matrix (n x k) modeled as fixed effects variables 
# @return NULL
# 
validate.vrPEER.data.objects <- function(Q, y, Z, X){
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(Q)[2] == dim(Z)[2])) stop("dim(L)[2] != dim(Z)[2]")
  if (!is.null(X)) if(!(dim(y)[1] == dim(X)[1])) stop("dim(y)[1] != dim(X)[1]")
}




#' Graph-constrained regression with variable-reduction procedure to handle the non-invertibility of 
#' a graph-originated penalty matrix
#' 
#' @param Q graph-originated penalty matrix \eqn{(p \times p)}, typically: a graph Laplacian matrix
#' @param y response values matrix \eqn{(n \times 1)}
#' @param Z design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling)
#' @param X design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling)
#' @param Q.normalized logical whether or not to use normalized version of a graph-originated penalty matrix
#' @param sv.thr - threshold value above which singular values of \code{Q} are considered "zeros" 
#' @param compute.Grace logical whether or not to compute \href{https://arxiv.org/abs/1506.08339}{Grace test} 
#' (a significance test for graph-constrained estimation) results 
#' 
#' @return 
#' \item{b.est}{vector with estimated values of \eqn{b} coefficients}
#' \item{beta.est}{vector with estimated values of \eqn{\beta} coefficients}
#' \item{lambda.Q}{value of \eqn{\lambda_Q} regularization parameter}
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
#' vrPEER.fit <- vrPEER(Q = L.norm, y = Y, Z = Z, X = NULL)
#'
#' # b coefficient estimates
#' vrPEER.fit$b.est
#' # regularization parameter values
#' vrPEER.fit$lambda.Q
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
#' vrPEER.fit <- vrPEER(Q = L.norm, y = Y, Z = Z, X = X)
#' 
# '# b coefficient estimates
#' vrPEER.fit$b.est
#' # beta coefficient estimates
#' vrPEER.fit$beta.est
#' # regularization parameters values
#' vrPEER.fit$lambda.Q
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
vrPEER <- function(Q, y, Z, X = NULL, 
                   Q.normalized = FALSE, 
                   sv.thr = 1e-5, 
                   compute.Grace = FALSE){

  # Check for data objects dimensions consistency 
  validate.vrPEER.data.objects(Q, y, Z, X)
  
  # Transform data objects 
  if (Q.normalized) Q <- L2L.normalized(Q)
  
  # Construct data objects used in variable-reduction approach
  pcov <- ifelse(is.null(X), 0, dim(X)[2])
  Z.ext <- cbind(X, Z) 
  Q.ext <- as.matrix(magic::adiag(matrix(0, nrow = pcov, ncol = pcov), Q))  # Q extension with 0-filled matrix to reflect penalty for predictors we do not want to penalize in estimation
  n <- dim(Z.ext)[1]
  p <- dim(Z.ext)[2]
  
  # SVD of Q 
  Q.ext.svd <- svd(Q.ext)
  U <- Q.ext.svd$u
  Sigma <- diag(Q.ext.svd$d)
  r <- length(which(diag(Sigma) < sv.thr))  # how many singular values of Sigma are "zeros"
  k <- p - r # how many singular values of Sigma are non-"zeros" (are OK)
  
  if (r == 0){
    # Define LMM random effects matrix lme.Z
    Q.chol <- chol(Q)
    Q.chol.inv <- solve(Q.chol)
    lme.Z <- Z %*% Q.chol.inv
    # Fit LMM
    id.bd1 <- factor(rep(1, n))
    lme.fit <- lme(fixed = y ~ 1, random = list(id.bd1 = pdIdent(~ lme.Z - 1)), method = "REML") 
  } else {
    # Procedure step (1): transform PEER to weighted Ridge
    Sigma.wave <- diag(diag(Sigma)[1:k])  # reduce Sigma to keep only "non-zeros" singular values 
    Z.wave <- Z.ext %*% U
    X.wave <- as.matrix(Z.wave[, (k+1):p])  # fixed effects matrix 
    Z.wave.red <- Z.wave[, 1:k]  # reduce Z.wave to keep only "non-zeros" singular values 
    # Procedure step (2): transform weighted Ridge to Ridge 
    b <- diag(sqrt(diag(Sigma.wave)))
    b.inv <- solve(b)
    Z.wave.red.wave <- Z.wave.red %*% b.inv
    # Fit LMM model
    id.bd1 <- factor(rep(1, n))
    lme.fit <- lme(fixed = y ~ X.wave + 1, random = list(id.bd1 = pdIdent(~ Z.wave.red.wave - 1)), method = "REML") 
  }
  
  # LMM-originated lambda parameter
  sigma.eps <- lme.fit$sigma 
  sigma.u <- as.numeric(as.matrix(VarCorr(lme.fit))[1, 2])
  lambda.Q <- (sigma.eps^2)/(sigma.u^2)
  
  # Compute Grace Test results
  if (!compute.Grace){
    grace.test.res <- NULL
  } else {
    grace.test.res <- tryCatch({
      Q.opt <- lambda.Q * Q.ext
      gt.res <- grace.test.C(Y = y, X = Z.ext, L = Q.opt, lambda.L = 1, lambda.2 = 0, normalize.L = F,
                             normalize.X  = FALSE, center.Y = FALSE)
      gt.beta <- gt.res$beta
      gt.beta.pval <- gt.res$pvalue
      if (pcov > 0){
        gt.b.est <- gt.beta[(pcov+1):(pcov+dim(Z)[2])]
        gt.b.est.pval <- gt.beta.pval[(pcov+1):(pcov+dim(Z)[2])]
      } else {
        gt.b.est <- gt.beta
        gt.b.est.pval <- gt.beta.pval
      }
      list(b.est = gt.b.est, b.est.pval = gt.b.est.pval, Q.opt = Q.opt)
    }, error = function(e) {
      msg <- "ERROR OCURRED WHEN PROCEEDING WITH grace.test"
      message(msg); message(e)
      msg
    }) 
  }
  
  # Compute coefficient estimates
  beta.b.est <- as.vector(solve(t(Z.ext) %*% Z.ext + lambda.Q * Q.ext) %*% t(Z.ext) %*% y)
  if (is.null(colnames(Z.ext))) names(beta.b.est) <- colnames(Z.ext)
  if (pcov > 0){
    b.est <- beta.b.est[(pcov+1):(pcov+dim(Z)[2])]
    beta.est <- beta.b.est[1:pcov]
  } else {
    b.est <- beta.b.est
    beta.est <- c()
  }
  
  res.list <- list(b.est = b.est, 
                   beta.est = beta.est, 
                   lambda.Q = lambda.Q,
                   grace.test.res = grace.test.res)
  return(res.list)
}