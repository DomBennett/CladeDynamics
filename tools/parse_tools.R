## 11/09/2014
## D.J. Bennett
## Tools for modifying trees

## Libraires
library (MoreTreeTools)

## Functions
safeChronos <- function (tree) {
  ## Wrapper for chronoMPL to handle unexpected errors
  ## see -- https://stat.ethz.ch/pipermail/r-sig-phylo/2014-April/003416.html
  temp <- try (chronoMPL (tree), silent = TRUE)
  if (any (class (temp) == 'phylo')) {
    tree <- temp
  }
  tree
}

convertToDist <- function  (tree, n = tree.dist) {
  # Convert a polytmous tree into a distribution of dichotomous trees
  trees <- list ()
  for (i in 1:n) {
    trees <- c (trees, list (multi2di (tree)))
  }
  class (trees) <- 'multiPhylo'
  trees
}