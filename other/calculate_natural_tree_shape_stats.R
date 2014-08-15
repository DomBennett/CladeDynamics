## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for natural trees

## Parameters
target <- 100
leeway <- 10

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Dirs
data.dir <- 'data'

## Input
natural.trees <- list ()
weights <- NULL # vector of source of tree, for weighting mean
tree.files <- list.files (data.dir, '\\.tre')
min.n <- target - (target*leeway/100)
max.n <- target + (target*leeway/100)
for (i in 1:length (tree.files)) {
  tree <- read.tree (file.path (data.dir, tree.files[i]))
  # choose the biggest tree if list
  if (class (tree) == 'multiPhylo') {
    sizes <- unlist(lapply (tree, getSize))
    tree <- tree[[which (sizes == max(sizes))[1]]]
  }
  # if the tree is bigger than target extract clades that are
  #  within leeway of target
  if (getSize (tree) >= max.n) {
    # may not work for polytomous trees, use try
    clade.trees <- try (getSubtrees (tree, min.n, max.n), silent = TRUE)
    if (is.null (clade.trees) | class (clade.trees) == 'try-error') {
      next
    }
    for (clade.tree in clade.trees) {
      natural.trees <- c (natural.trees, list (clade.tree))
      weights <- c (weights, i)
    }
  } else if (getSize (tree) >= min.n) {
    natural.trees <- c (natural.trees, list (tree))
    weights <- c (weights, i)
  } else {
    next
  }
}

## Calculate tree stats for each group of trees
colless.stat <- sackin.stat <- iprime.stat <- gamma.stat <- tc.stat <- NULL
for (w in unique (weights)) {
  natural.tree.stats <- calcTreeShapeStats (natural.trees[weights == w])
  colless.stat <- c (colless.stat, natural.tree.stats['mean.colless.stat'][[1]])
  sackin.stat <- c (sackin.stat, natural.tree.stats['mean.sackin.stat'][[1]])
  iprime.stat <- c (iprime.stat, natural.tree.stats['mean.iprime.stat'][[1]])
  gamma.stat <- c (gamma.stat, natural.tree.stats['mean.gamma.stat'][[1]])
  tc.stat <- c (tc.stat, natural.tree.stats['mean.tc.stat'][[1]])
}
natural.tree.stats <- data.frame (colless.stat, sackin.stat, iprime.stat,
                                  gamma.stat, tc.stat)

## Save tree stats
filename <- paste0 ('natural_tree_stats_t', target, '_l', leeway, '.Rd')
save (natural.tree.stats, file = file.path (data.dir, filename))