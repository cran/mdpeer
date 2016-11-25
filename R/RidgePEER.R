
# Comoute optimization problem objective function value.
# 
# @param lambda function parameter
# @param S function parameter
# @param y function parameter
# @param Z function parameter
# 
# @return objective function value
# 
objective.fun.h <- function(lambda, S, y, Z){
  lambda.Q <- lambda[1]
  lambda.R <- lambda[2]
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  D <- lambda.Q * S + lambda.R * diag(p)
  val <- n * log(sum(y^2) -  t(y) %*% Z %*% solve(D + t(Z) %*% Z) %*% t(Z) %*% y) + log(det((D + t(Z) %*% Z) %*% solve(D)))
  return(val)
}



# Wrapper to allow for adding additional parameters to h and gradient of h functions.
# 
# @param func function to be wrapped
# 
# @return wrapped function
# 
applyDefaults <- function(func, ...){
  function(x) func(x, ...)
}



# Perform custom optimization with the use of one of the optim function algorithms provided.
# 
# @param optim.function optimization function signature
# @param x0 parameter initial values
# @param fn function to be minimized
# @param lower lower boundary for parameter values
# @param upper upper boundary for parameter values
# 
# @return optimizaton function result object
# 
optimize.custom <- function(optim.function, x0, fn, lower, upper){
  
  if (optim.function == "sbplx"){
    opt.out <- sbplx(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
    
  } else if (optim.function == "bobyqa"){
    opt.out <- bobyqa(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
    
  } else if (optim.function == "cobyla"){
    opt.out <- cobyla(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
    
  } else if (optim.function == "lbfgs"){
    opt.out <- lbfgs(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
    
  } else if (optim.function == "neldermead"){
    opt.out <- neldermead(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
    
  } else if (optim.function == "auglag"){
    opt.out <- auglag(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
    
  } else if (optim.function == "slsqp"){
    opt.out <- slsqp(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
    
  } else if (optim.function == "varmetric"){
    opt.out <- varmetric(x0, fn, lower = lower, upper = upper) 
    opt.lambda.Q <- opt.out$par[1]
    opt.lambda.R <- opt.out$par[2]
    opt.bestval <- opt.out$value
  }
  
  opt.res.list <- list(opt.out = opt.out, opt.bestval = opt.bestval,
                       opt.lambda.Q = opt.lambda.Q, opt.lambda.R = opt.lambda.R)
  return(opt.res.list)
}




#' @title Graph-constrained estimation with enhanced regulariazation parameter selection 
#' 
#' @description 
#' Perform penalized estimation where penalty term is a linear combination of 
#' graph-originated penalty and Ridge-originated penalty terms. Corresponding regularization
#' parameter values \eqn{\lambda} (\eqn{\lambda_Q, \lambda_R}) are estimated as ML 
#' estimators from corresponding Linear Mixed Model solution.
#' 
#' @details
#' We assume a model: \deqn{y = X\beta + Zb + \varepsilon} and we 
#' estimate its coefficients  \eqn{\beta, b} as follows: 
#' \deqn{\widehat{\beta}, \widehat{b}= arg \; min_{\beta,b}\left \{ (y - X\beta - Zb)^T(y - X\beta - Zb) + \lambda_Qb^TQb +  \lambda_Rb^Tb\right \}, }
#' or, equivalently: 
#' \deqn{\widehat{\beta}, \widehat{b}= arg \; min_{\beta,b}\left \{ (y - X\beta - Zb)^T(y - X\beta - Zb) + \lambda_Q (b^TQb +  \lambda_2 b^Tb\right) \}, }
#'   where: 
#' \itemize{
#'  \item \eqn{y} - a response variable,  
#'  \item \eqn{X} - matrix of variables whose coefficients we do not want to penalize in the estimation process (assumed to be standarized to have mean equal 0 and variance equal 1),  
#'  \item \eqn{Z} - matrix of variables whose coefficients we want to penalize in the estimation process (assumed to be standarized to have mean equal 0 and variance equal 1),  
#'  \item \eqn{Q} - a graph-originated penalty matrix, typically: a graph Laplacian matrix, 
#'  \item \eqn{\lambda_Q, \lambda_R, \lambda_2} - regularization parameters (\eqn{\lambda_R = \lambda_Q  \lambda_2}). 
#'} 
#' 
#' In the coefficient estimation formula, \eqn{b^TQb} denotes a graph-originated model penalty term, 
#' \eqn{b^Tb} denotes a Ridge-originated model penalty term.
#' 
#' A graph-originated penalty term allows imposing similarity between coefficients 
#' of variables which are similar (or connected), based on some graph information given. 
#' Additional Ridge-originated penalty term facilitates parameters estimation: 
#' it reduces computational issues (arising from singularity in a graph-originated penalty matrix) 
#' and yields plausible results in situations when graph information is not informative or when 
#' it is unclear whether connectivities represented by a graph reflect similarities among corresponding coefficients.
#' 
#' Implementation of model regularization parameters \eqn{\lambda} (\eqn{\lambda_Q, \lambda_R}) selection 
#' utilizes the known fact of equivalence between penalized regression and Linear Mixed Model solutions, 
#' and provides values of the two regularization parameters that are Maximum Likelihood estimators of the latter.
#' 
#' Several optimization algorithms for minimization objective function in MLE estimation 
#' are available: \code{sbplx}, \code{bobyqa}, \code{cobyla}, \code{lbfgs}, \code{neldermead}, 
#' \code{auglag}, \code{slsqp}, \code{varmetric}. Simulation study conducted shown that \code{sbplx} 
#' algorithm was of the best performance in the cases investigated. This algorithm is set as default
#' optimization algorithm choice and is strongly recommended by the package authors for 
#' the \code{RidgePEER} function. All the optimization algorithms available are implemented
#' in the \code{nloptr} package. Please refer to the \code{nloptr} package for 
#' details. 
#' 
#' @param Q graph-originated penalty matrix \eqn{(p \times p)}, typically: a graph Laplacian matrix
#' @param y response values matrix \eqn{(n \times 1)}
#' @param Z design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling)
#' @param X design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling);
#' note that an additional column representing intercept will always be added to \code{X} for the sake of computational simplicity 
#' (particularly, if \code{X} is initially \code{NULL}, it will be transformed so as to contain a single column of 1's)
#' @param Q.normalized logical whether or not to use normalized version of a graph-originated penalty matrix
#' @param add.Ridge logical whether or not to include Ridge penalty term in \code{RidgePEER} model penalty term 
#' @param add.PEER logical whether or not to include PEER penalty term in \code{RidgePEER} model penalty term 
#' @param optim.function signature of optimization algorithm used in MLE estimation (see: Details)
#' @param x0 2-elements vector of initial values for \eqn{\lambda_Q, \lambda_R} parameters MLE 
#' estimation in optimization algorithm 
#' @param lambda.Q.lo lower boundary for values space in which we search for optimal value of \eqn{\lambda_Q}
#' @param lambda.Q.up upper boundary for values space in which we search for optimal value of \eqn{\lambda_Q}
#' @param lambda.R.lo lower boundary for values space in which we search for optimal value of \eqn{\lambda_R}
#' @param lambda.R.up upper boundary for values space in which we search for optimal value of \eqn{\lambda_R}
#' @param verbose logical whether or not message out information from function execution
#' @param compute.Grace logical whether or not to compute \href{https://arxiv.org/abs/1506.08339}{Grace test} 
#' (a significance test for graph-constrained estimation) results 
#'
#' @return 
#' \item{b.est}{vector with estimated values of \eqn{b} coefficients}
#' \item{beta.est}{vector with estimated values of \eqn{\beta} coefficients}
#' \item{opt.out}{output from the \eqn{\lambda_Q, \lambda_R} parameters optimization algorithm}
#' \item{opt.bestval}{objective function value from the \eqn{\lambda_Q, \lambda_R} parameters optimization algorithm}
#' \item{lambda.Q}{ML estimator value of \eqn{\lambda_Q} regularization parameter}
#' \item{lambda.R}{ML estimator value of \eqn{\lambda_2} regularization parameter}
#' \item{lambda.2}{\code{lambda.R}/\code{lambda.Q} value}
#' \item{grace.test.res}{output from the Grace significance test for graph-constrained estimation computation}
#' 
#' @examples 
#' # Example 1. model without covariates
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
#' RidgePEER.fit <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = NULL)
#' 
#' # b coefficient estimates
#' RidgePEER.fit$b.est
#' # MLE of regularization parameters lambda
#' RidgePEER.fit$lambda.Q
#' RidgePEER.fit$lambda.R
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
#' RidgePEER.fit <- RidgePEER(Q = L.norm, y = Y, Z = Z, X = X)
#' 
#' # b coefficient estimates
#' RidgePEER.fit$b.est
#' # beta coefficient estimates
#' RidgePEER.fit$beta.est
#' # MLE of regularization parameters lambda
#' RidgePEER.fit$lambda.Q
#' RidgePEER.fit$lambda.R
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
#' 
#'
#' @import nloptr
#' @export
#' 
RidgePEER <- function(Q, y, Z, X = NULL,
                      Q.normalized = FALSE,
                      add.Ridge = TRUE,
                      add.PEER = TRUE,
                      optim.function = "sbplx",
                      x0 = c(1, 1),
                      lambda.Q.lo = 10^(-5),
                      lambda.Q.up = 1e6,
                      lambda.R.lo = 10^(-5),
                      lambda.R.up = 1e6,
                      verbose = TRUE,
                      compute.Grace = FALSE){
  
  # Check for data objects dimensions consistency 
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(Q)[2] == dim(Z)[2])) stop("dim(Q)[2] != dim(Z)[2]")
  if (!is.null(X)) if(!(dim(y)[1] == dim(X)[1])) stop("dim(y)[1] != dim(X)[1]")
  if ((!add.Ridge) & (!add.PEER)) stop("(!add.Ridge) & (!add.PEER)")
  
  # Transform data objects / parameters 
  if (Q.normalized) Q <- L2L.normalized(Q)
  if (!add.Ridge){
    lambda.R.lo <- 10^(-5)
    lambda.R.up <- 10^(-5)*2
    x0[2] <- 10^(-5)*1.5
  }
  if (!add.PEER){
    Q <- diag(dim(Z)[2])
  }
  
  # Define data object dimensions
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  
  # Add intercept to X design matrix 
  if (!is.null(X)){
    X <- cbind(1, X)
  } else {
    X <- matrix(rep(1, n), ncol = 1)
  }
  
  # Regress out fixed effect coefficients 
  X.hash <- diag(dim(X)[1]) - X %*% solve(t(X) %*% X) %*% t(X)
  y.wave <- X.hash %*% y
  Z.wave <- X.hash %*% Z
  
  # Transform Laplacian-constrained problem into generalized Ridge-constrained problem 
  Q.svd <- svd(Q)
  S <- diag(Q.svd$d)
  V <- Q.svd$u
  Z.waveV <- Z.wave %*% V

  
  # ------------------------------------------------------------------------------
  
  # Define objective function for data objects 
  fn <- applyDefaults(objective.fun.h, S = S, y = y.wave, Z = Z.waveV)
  
  # Define lower and upper search boundaries
  lower <- c(lambda.Q.lo, lambda.R.lo)
  upper <- c(lambda.Q.up, lambda.R.up) 
  
  # Optimize 
  if(verbose) message(paste0("Running optimizer: ", optim.function, "..."))
  opt.res.list <- optimize.custom(optim.function, x0, fn, lower, upper)
  opt.out <- opt.res.list$opt.out
  opt.bestval <- opt.res.list$opt.bestval
  opt.lambda.Q = opt.res.list$opt.lambda.Q
  opt.lambda.R = opt.res.list$opt.lambda.R
  if(verbose) message(paste0("Model lambdas optimized. Optimal values:",
                             "\nopt.lambda.Q: ", round(opt.lambda.Q, 3), 
                             "\nopt.lambda.R: ", round(opt.lambda.R, 3)))
  
  # Compute Grace Test results
  if (!compute.Grace){
    grace.test.res <- NULL
  } else {
    grace.test.res <- tryCatch({
      Q.opt <- opt.lambda.Q * Q + opt.lambda.R * diag(p)
      gt.res <- grace.test.C(Y = y.wave, X = Z.wave, L = Q.opt, lambda.L = 1, lambda.2 = 0, normalize.L = F,
                             normalize.X  = FALSE, center.Y = FALSE)
      list(b.est = gt.res$beta, b.est.pval = as.vector(gt.res$pvalue), Q.opt = Q.opt)
    }, error = function(e) {
      msg <- "ERROR OCURRED WHEN PROCEEDING WITH grace.test"
      message(msg); message(e)
      msg
    }) 
  }
  
  # Define etimated b, beta coefficient vectors 
  ## Define penalty term D with the use of optimal lambdas values
  opt.D <- opt.lambda.Q * S + opt.lambda.R * diag(p)
  ## Estimate random effect coefficients from assumed generalized Ridge-constrained problem
  b.est.gR <- as.vector(solve(t(Z.waveV) %*% Z.waveV + 1 * opt.D) %*% t(Z.waveV) %*% y.wave)
  ## Convert estimates from generalized Ridge-constrained problem into Laplacian-constrained problem 
  b.est <- V %*% b.est.gR
  ## Estimate fixed effect coefficients 
  beta.est <- as.vector(solve(t(X) %*% X) %*% t(X) %*% (y - Z %*% b.est))
  names(beta.est)[1] <- "Intercept"
  # if(!is.null(colnames(X)) & dim(X)[2]>1) (beta.est)[2:dim(X)[2]] <- colnames(X)
  # if(!is.null(colnames(Z))) names(b.est) <- colnames(Z)
  
  res.list <- list(b.est = b.est, 
                   beta.est = beta.est,
                   opt.out = opt.out, 
                   opt.bestval = opt.bestval,
                   lambda.Q = opt.lambda.Q,
                   lambda.R = opt.lambda.R,
                   lambda.2 = opt.lambda.R / opt.lambda.Q,
                   grace.test.res = grace.test.res)
  return(res.list)
}

