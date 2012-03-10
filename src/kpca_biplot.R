library(svd)
library(kernlab)

# Take a matrix X and a parameter a and produce a list of the matrices G and H
# that form its singular value decomposition as done in the paper.
our.svd <- function(X, a) {
  decomp <- svd(X)
  U <- decomp$u
  V <- decomp$v
  D <- diag(decomp$d)
  return(list(G = U %*% D^a, H = V %*% D^(1-a)))
}

# Compute the kpca biplot, with the kernlab kpca package doing most of the work.
# This function is not yet finished. I still need to add in the biplot part of
# kpca biplot.
kpca.biplot <- function(x, kernel = "rbfdot", kpar = list(sigma = 0.1),
                        features = 0, th = 1e-4, na.action = na.omit,
                        alpha = 0.5, ...)
  {
    x.svd <- our.svd(x, alpha)
    G <- x.svd$G
    H <- x.svd$H

    res <- kpca(t(H), kernel, kpar, features, th, na.action)
    kernel <- do.call(kernel, kpar)
    kg <- kernelMatrix(kernel, G, t(H))

    # Projects G onto the principal component values
    G.prj <- kg %*% pcv(res)

    # Right now it just returns a list that contains the rows of H and the rows
    # of G projected onto the principal components.
    return(list(g = G.prj, h = rotated(res)))
  }
    
