## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for natural trees

## Parameters
max.n = 1000
min.n = 500

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Functions
getBigClades <- function (tree, min.n, max.n) {
  countChildren <- function (node) {
    length (getChildren (tree, node))
  }
  ntips <- getSize (tree)
  nodes <- ntips:(ntips + tree$Nnode)
  nchildren <- mdply (.data = data.frame (node = nodes),
                      .fun = countChildren)
  bool <- nchildren$V1 > min.n & nchildren$V1 < max.n
  nchildren$node[bool]
}

## Dirs
data.dir <- 'data'

## Input
natural.trees <- list ()
tree.files <- list.files (data.dir, '.tre')
for (i in 1:length (tree.files)) {
  tree <- read.tree (file.path (data.dir, tree.files[i]))
  # choose the first tree if list
  if (class (tree) == 'multiPhylo') {
    tree <- tree[[1]]
  }
  # if the tree is bigger than 1000 extract clades that are
  #  greater than 500 and less than 1000
  if (getSize (tree) > 1000) {
    nodes <- getBigClades (tree, min.n, max.n)
    for (n in nodes) {
      clade.tree <- extract.clade (tree, n)
      natural.trees <- c (natural.trees, list (clade.tree))
    }
  } else {
    natural.trees <- c (natural.trees, list (tree))
  }
}

## Calculate tree stats
natural.tree.stats <- calcTreeShapeStats (natural.trees)

## Save tree stats
save (natural.tree.stats, file = file.path (data.dir, 'natural_tree_stats.Rd'))