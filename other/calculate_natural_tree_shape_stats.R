## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for natural trees

## Parameters
target <- 100
leeway <- 10
tree.dist <- 1 # the number of trees in a distribution for a polytomous tree

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Functions
safeChronos <- function (tree, ...) {
  ## Wrapper for chronos to handle unexpected errors
  ## see -- https://stat.ethz.ch/pipermail/r-sig-phylo/2014-April/003416.html
  temp <- try (chronos (tree, ...), silent = TRUE)
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
package <- function (tree, max.n, min.n) {
  # Take tree or list of trees, return phylo or multiphylo of right sized trees
  .resize <- function (tree) {
    if (getSize (tree) >= max.n) {
      # if the tree is bigger than target extract clades that are
      #  within leeway of target
      clade.trees <- getSubtrees (tree, min.n, max.n)
      if (!is.null (clade.trees)) {
        if (any (class (clade.trees) == 'multiPhylo')) {
          return (list (clade.trees))
        } else {
          return (list (clade.trees))
        }
      }
    } else if (getSize (tree) >= min.n) {
      # if not... just add it
      return (list (tree))
    }
  }
  trees <- list ()
  # tree is a multiPhylo, run .add for each
  if (any (class (tree) == 'multiPhylo')) {
    clade.trees <- list ()
    for (t in tree) {
      temp <- .resize (t)
      if (!is.null (temp)) {
        clade.trees <- c (clade.trees, temp)
      }
    }
    if (length (clade.trees) > 0) {
      trees <- c (trees, list (clade.trees))
    }
  } else {
    temp <- .resize (tree)
    if (!is.null (temp)) {
      return (temp)
    }
  }
  if (length (trees) > 1) {
    class (trees) <= 'multiPhylo'
    return (trees)
  }
  if (length (trees) == 1) {
    trees <- trees[[1]]
    return (trees)
  }
}

## Dirs
data.dir <- file.path ('data', 'trees')

## Input
trees <- list () # hold sets of natural trees
tree.files <- list.files (data.dir, '\\.tre') # list all .tre files in data folder
tree.files <- sample (tree.files, 10)
min.n <- target - (target*leeway/100) # smallest ...
max.n <- target + (target*leeway/100) # ... and biggest sizes for trees
for (i in 1:length (tree.files)) {
  # suppress add.terminal warnings
  tree <- suppressWarnings (read.tree (file.path (data.dir, tree.files[i])))
  # choose the biggest tree of the trees for a study
  if (class (tree) == 'multiPhylo') {
    sizes <- unlist(lapply (tree, getSize))
    tree <- tree[[which (sizes == max(sizes))[1]]]
  }
  if (getSize (tree) < min.n) {
    next
  }
  # does it have branch lengths?
  bl.bool <- !is.null (tree$edge.length) && all (!is.na (tree$edge.length))
  if (!bl.bool) {
    # make sure if any edgelenths are NA, edgelenghts are removed
    tree$edge.length <- NULL
  }
  # is it ultrametric?
  ultra.bool <- bl.bool && is.ultrametric (tree)
  # is it polytomous?
  poly.bool <- getSize (tree) != (tree$Nnode + 1)
  # print progress
  cat (paste0 ('\nWorking on [', tree.files[i],']'))
  cat (paste0 ('\n.... [', getSize (tree), '] taxa'))
  cat (paste0 ('\n.... [', bl.bool, '] branch lengths'))
  cat (paste0 ('\n.... [', ultra.bool, '] ultrametric'))
  cat (paste0 ('\n.... [', !poly.bool, '] dichtomous'))
  # if not ultrametric make it (if I can)
  if (bl.bool && !ultra.bool) {
    tree <- safeChronos (tree, lambda = 1, quiet = TRUE)
  }
  # if polytomous, convert to a distribution
  if (poly.bool) {
    tree <- convertToDist (tree)
  }
  # add tree to list of trees
  res <- package (tree, max.n, min.n)
  if (!is.null (res)) {
    trees <- c (trees, res)
  }
}

## Calculate tree stats for each set of trees
colless.stat <- sackin.stat <- iprime.stat <- gamma.stat <- tc.stat <- NULL
for (set in trees) {
  ## TODO -- fix the double class problem: chronos and phylo
  stats <- calcTreeShapeStats (set) # calculate all stats
  # extract the mean value of the set
  colless.stat <- c (colless.stat, stats['mean.colless.stat'][[1]])
  sackin.stat <- c (sackin.stat, stats['mean.sackin.stat'][[1]])
  iprime.stat <- c (iprime.stat, stats['mean.iprime.stat'][[1]])
  gamma.stat <- c (gamma.stat, stats['mean.gamma.stat'][[1]])
  tc.stat <- c (tc.stat, stats['mean.tc.stat'][[1]])
}
natural.tree.stats <- data.frame (colless.stat, sackin.stat, iprime.stat,
                                  gamma.stat, tc.stat)

## Save tree stats
filename <- paste0 ('natural_tree_stats_t', target, '_l', leeway, '.Rd')
save (natural.tree.stats, file = file.path (data.dir, filename))