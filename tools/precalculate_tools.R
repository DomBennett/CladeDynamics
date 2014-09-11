## 11/09/2014
## D.J. Bennett
## Tools for precalculating tree shape stats

## Libraries
library (MoreTreeTools)

## Functions
pack <- function (treeobj, max.n, min.n) {
  # Take phylo or multiphylo, return multiphylo of 
  #  right-sized trees
  if (any (class (treeobj) == 'multiPhylo')) {
    trees <- list ()
    class (trees) <- 'multiPhylo'
    for (tree in treeobj) {
      tree <- getSubtrees (tree, min.n = min.n,
                           max.n = max.n)
      if (!is.null (tree)) {
        trees <- cTrees (trees, tree)
      }
    }
    if (length (trees) > 0) {
      return (trees)
    } else {
      # here no subtrees have been extracted
      return (NULL)
    }
  } else if (any (class (treeobj) == 'phylo')) {
    tree <- getSubtrees (treeobj, min.n = min.n,
                         max.n = max.n)
    return (tree)
  } else {
    cat ('Unknown tree object given: not phylo or multiphylo')
    return (NULL)
  }
}