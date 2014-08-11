## 11/08/2014
## D.J. Bennett
## Script for downloading species trees from treebase

## Libraries
library (treebase)

## Search treebase for data on all trees
meta <- metadata ()

## Download all species trees with more than 100 taxa
downloaded.trees <- list ()
study.data <- list ()
for (i in 1:nrow (meta)) {
  tree.data <- meta[i, ]
  if (tree.data$kind == 'Species Tree') {
    if (tree.data$ntaxa >= 100) {
      tree <- search_treebase(tree.data$Tree.id, by = 'id.tree')
      if (length (tree) > 0) {
        downloaded.trees <- c (downloaded.trees, tree)
        study.data <- c (study.data, list (tree.data))
      }
    }
  }
}
