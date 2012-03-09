library(svd)

# Take a matrix X and a parameter a and produce a list of the matrices G and H
# that form its singular value decomposition as done in the paper.
our.svd <- function(X, a) {
  decomp <- svd(X)
  U <- decomp$u
  V <- decomp$v
  D <- diag(decomp$d)
  return(list(G = U %*% D^a, H = V %*% D^(1-a)))
}
