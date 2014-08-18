## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for natural trees

## Parameters
target <- 100
leeway <- 10

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Functions
convertToDist <- function  (tree, n = 100) {
  # Convert a polytmous tree into a distribution of dichotomous trees
  trees <- list ()
  for (i in 1:n) {
    trees <- c (trees, list (multi2di (tree)))
  }
  class (trees) <- 'multiPhylo'
  trees
}
addToList <- function (tree, trees, max.n, min.n) {
  # Take tree and list, return list with right sized tree added
  .resize <- function (tree) {
    if (getSize (tree) >= max.n) {
      # if the tree is bigger than target extract clades that are
      #  within leeway of target
      clade.trees <- getSubtrees (tree, min.n, max.n)
      if (!is.null (clade.trees)) {
        if (class (clade.trees) == 'multiPhylo') {
          return (clade.trees)
        } else {
          return (list (clade.trees))
        }
      }
    } else if (getSize (tree) >= min.n) {
      # if not... just add it
      return (list (tree))
    }
  }
  # tree is a multiPhylo, run .add for each
  if (class (tree) == 'multiPhylo') {
    clade.trees <- list ()
    for (t in tree) {
      clade.trees <- c (clade.trees, .resize (t))
    }
    trees <- c (trees, list (clade.trees))
  } else {
    trees <- c (trees, .resize (tree))
  }
  trees
}

## Dirs
data.dir <- 'data'

## Input
trees <- list () # hold sets of natural trees
tree.files <- list.files (data.dir, '\\.tre') # list all .tre files in data folder
min.n <- target - (target*leeway/100) # smallest ...
max.n <- target + (target*leeway/100) # ... and biggest sizes for trees
for (i in 1:length (tree.files)) {
  tree <- read.tree (file.path (data.dir, tree.files[i]))
  # choose the biggest tree of the trees for a study
  if (class (tree) == 'multiPhylo') {
    sizes <- unlist(lapply (tree, getSize))
    tree <- tree[[which (sizes == max(sizes))[1]]]
  }
  # if polytomous, convert to a distribution
  if (getSize (tree) != (tree$Nnode + 1)) {
    tree <- convertToDist (tree)
  }
  # add tree to list of trees
  trees <- addToList (tree, trees, max.n, min.n)
}

## Calculate tree stats for each set of trees
colless.stat <- sackin.stat <- iprime.stat <- gamma.stat <- tc.stat <- NULL
for (set in trees) {
  if (length (set) == 0) {
    # if no trees met criteria, then set will be empty
    next
  }
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