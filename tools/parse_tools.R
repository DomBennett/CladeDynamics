## 11/09/2014
## D.J. Bennett
## Tools for modifying trees

## Libraires
library (MoreTreeTools)

## Functions
runRateSmoother <- function (trees, rsmoother, i) {
  if (rsmoother == 'pathD8') {
    rsmoother <- pathD8
  } else if (rsmoother == 'chronoMPL') {
    rsmoother <- safeChronoMPL
  } else if (rsmoother == 'chronopl') {
    rsmoother <- safeChronopl
  } else if (rsmoother == 'chronos') {
    rsmoother <- safeChronos
  } else {
    stop ('Rate smoother must be pathD8, chronoMPL, chronos or chronopl')
  }
  # run rate smoother over multiple trees
  .run <- function (j) {
    tree <- rsmoother (trees[[j]], i)
    new.trees <<- c (new.trees, list (tree))
  }
  if (class (trees) == 'multiPhylo') {
    new.trees <- list ()
    m_ply (.data=data.frame(j=1:length(trees)), .fun=.run)
    class (new.trees) <- 'multiPhylo'
    return (new.trees)
  } else {
    return (rsmoother (trees, i))
  }
}

safeChronoMPL <- function (tree, i=1) {
  ## Wrapper for chronoMPL to handle unexpected errors
  ## see -- https://stat.ethz.ch/pipermail/r-sig-phylo/2014-April/003416.html
  temp <- try (chronoMPL (tree), silent = TRUE)
  if (any (class (temp) == 'phylo')) {
    tree <- temp
  }
  tree
}

safeChronos <- function (tree, i=1) {
  if (!is.ultrametric (tree)) {
    temp <- try (chronos (tree, quiet=TRUE), silent = TRUE)
    if (any (class (temp) == 'phylo')) {
      tree <- temp
    }
  }
  tree
}

safeChronopl <- function (tree, i=1) {
  temp <- try (chronopl (tree, lambda=1), silent = TRUE)
  if (any (class (temp) == 'phylo')) {
    tree <- temp
  }
  tree
}

pathD8 <- function (tree, i=1) {
  if (!is.ultrametric (tree)) {
    # Run pathd8 from system path
    # Use i to give unique name (optional)
    tree$node.label <- NULL  # reduce error rate by removing node labels that may be non-regular
    infile <- paste0 ('temp_pathd8_', i, '_in.tre')
    outfile <- paste0 ('temp_pathd8_', i, '_out.tre')
    write.tree (tree, infile)
    system (paste0 ('./PATHd8 ', infile, ' ', outfile),
            ignore.stdout=TRUE)
    tree <- read.tree (outfile)
    system (paste0 ('rm ', infile, ' ', outfile))
    tree <- tree[[1]]  # use d8 tree
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