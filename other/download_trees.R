## 11/08/2014
## D.J. Bennett
## Script for downloading species trees from treebase

## Libraries
library (treebase)

## Search treebase for data on all trees
meta <- metadata ()

## First clear out all treebase trees from dir
treebase.files <- list.files (path = 'data', pattern = 'treebase')
file.remove(file.path ('data', treebase.files))

## Download all species trees with more than 100 taxa
## Write the trees to file using their study_id
for (i in 1:nrow (meta)) {
  tree.data <- meta[i, ]
  if (tree.data$kind == 'Species Tree') {
    if (tree.data$ntaxa >= 100) {
      tree <- search_treebase(tree.data$Tree.id, by = 'id.tree')
      if (length (tree) > 0) {
        filename <- paste0 ('treebase_', tree.data$Study.id, '.tre')
        write.tree (tree[[1]], file = file.path ('data', filename),
                    append = TRUE)
      }
    }
  }
}