## 11/08/2014
## D.J. Bennett
## Stage: download real trees from TreeBase
## Note. There are lots of warnings, this is
##  not me but TreeBase

## Libraries
library (treebase)

## Dirs
output.dir <- file.path ('data', 'trees', 'treebase')

## Process
# First clear out all treebase trees from dir
treebase.files <- list.files (path = output.dir,
                              pattern = 'treebase')
file.remove(file.path (output.dir, treebase.files))
# setup download log
download.log <- file.path (output.dir, 'download_log.csv')
if (file.exists (download.log)) {
  file.remove (download.log)
}
headers <- data.frame ("Study.id", "Tree.id", "kind", "type",
                       "quality", "ntaxa", "date", "publisher",
                       "author", "title" )
write.table (headers, download.log, sep = ',', row.names = FALSE,
             col.names = FALSE)
# Search treebase for data on all trees
meta <- metadata ()
# Download all species trees with more than 100 taxa
# Write the trees to file using their study_id and log their data
for (i in 1:nrow (meta)) {
  tree.data <- meta[i, ]
  if (tree.data$kind == 'Species Tree') {
    if (tree.data$ntaxa >= 100) {
      tree <- search_treebase(tree.data$Tree.id, by = 'id.tree')
      if (length (tree) > 0) {
        filename <- paste0 ('treebase_', tree.data$Study.id, '.tre')
        write.tree (tree[[1]], file = file.path ('data', filename),
                    append = TRUE)
        write.table (tree.data, download.log, sep = ',', append = TRUE,
                     col.names = FALSE, row.names = FALSE)
      }
    }
  }
}