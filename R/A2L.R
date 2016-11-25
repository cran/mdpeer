#' Compute a graph Laplacian matrix from a graph Adjacency matrix 
#' 
#' @param Adj graph Adjacency matrix (squared symmetric matrix)
#' @return graph Laplacian matrix 
#' @examples 
#' # Create exemplary graph Adjacency matrix
#' p1 <- 10
#' p2 <- 90
#' p <- p1 + p2
#' A <- matrix(rep(0, p * p), p, p)
#' A[1:p1, 1:p1] <- 1
#' A[(p1 + 1):p, (p1 + 1):p] <- 1
#' vizu.mat(A, "Adjacency matrix")
#' 
#' # Compute corresponding graph Laplacian matrix
#' L <- Adj2Lap(A)
#' vizu.mat(L, "Laplacian matrix")
#' @export
#' 
Adj2Lap <- function(Adj){
  # D degree matrix: the sum of the weights of the connectivity matrix
  D.diag <- apply(Adj, MARGIN = 1, function(row) sum(row))
  D <- diag(D.diag, nrow = nrow(Adj), ncol = ncol(Adj))
  # Laplacian matrix 
  L <- D - Adj
  return(L)
}