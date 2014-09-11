## 23/07/2014
## D.J. Bennett
## Tools for testing

# Deps
library (plyr)
library (ape)

# Functions
genRandomData <- function (n.clades, time.span, parent.lambda) {
  # Generate random clade success data using Poisson dist
  .each <- function (i) {
    random.lambda <- rpois (1, 1)
    c (sort (rpois (n, random.lambda)),
       sort (rpois (n, random.lambda), decreasing = TRUE))
  }
  n <- round (time.span/2)
  res <- mdply (.data = data.frame (i = 1:n.clades), .fun = .each)[ ,-1]
  res <- t (res)
  colnames (res) <- paste0 ('n', 1:n.clades)
  rownames (res) <- 1:time.span
  res
}

genRandomTrees <- function (time.span) {
  # Generate a series of left growing trees
  trees <- list ()
  .each <- function (i) {
    tree <- stree (i, 'left')
    tree <- compute.brlen (tree)
    tree$node.label <- paste0 ('n', 1:tree$Nnode)
    trees <<- c (trees, list (tree))
  }
  m_ply (.data = data.frame (i = 3:(time.span + 2)), .fun = .each)
  trees
}