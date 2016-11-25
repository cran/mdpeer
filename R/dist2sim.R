#' Transform Distance matrix to Similarity matrix 
#' 
#' Transform Distance matrix into Similarity matrix. In particular, might be used
#' to obtain a graph Adjacency matrix from a graph Distance matrix. 
#' 
#' @param dist.mat distance matrix (squared symmetric matrix)
#' 
#' @return  similarity matrix 
#' 
#' @examples 
#' set.seed(123)
#' x <- matrix(rnorm(100), nrow = 10)
#' dist.mat <- as.matrix(dist(x))
#' vizu.mat(dist.mat, "Distance matrix")
#' 
#' sim.mat <- dist2sim(dist.mat)
#' vizu.mat(sim.mat, "Similarity matrix")
#' 
#' @export
#' 
dist2sim <- function(dist.mat){
  sim.mat <- 1/(dist.mat+1)
  diag(sim.mat) <- rep(0, dim(sim.mat)[2])
  return(sim.mat)
}