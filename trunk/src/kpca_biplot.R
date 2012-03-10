library(kernlab)
library(svd)

# Take a matrix X and a parameter a and produce a list of the matrices G and H
# that form its singular value decomposition G %*% t(H) as done in the paper.
our.svd <- function(X, a) {
  decomp <- svd(X)
  U <- decomp$u
  V <- decomp$v
  D <- diag(decomp$d)
  return(list(G = U %*% D^a, H = V %*% D^(1-a)))
}

# Computes the svd of x, breaking it into the components G and H. Does kpca
# on H and returns a list of G and H projected onto the principal
# components.
our.kpca <- function(x, kernel = "rbfdot", kpar = list(sigma = 0.1),
                        features = 0, th = 1e-4, na.action = na.omit,
                        alpha = 0.5, ...)
  {
    # Do all of the preprocessing nonsense.
    x <- as.matrix(x)
    x <- log2(x)
    x <- t(scale(t(x)))
    mn <- colMeans(x)
    x <- sweep(x, 2, mn)
    
    x.svd <- our.svd(x, alpha)
    G <- x.svd$G
    H <- x.svd$H

    res <- kpca(H, kernel, kpar, features, th, na.action)
    kernel <- do.call(kernel, kpar)
    kg <- kernelMatrix(kernel, H, G)
  
    # Projects G onto the principal component values
    G.prj <- t(kg) %*% pcv(res)

    # Right now it just returns a list that contains the rows of H and the rows
    # of G projected onto the principal components.
    return(list(micro.arrays = G.prj, gene.expressions = rotated(res)))
  }
    
