## 11/08/2014
## D.J. Bennett
## Stage: download real trees from TreeBase
## Note. There are lots of warnings, from
##  treebase functions, I've supressed them
## Note. 'Error: 1: Extra content at the end of the document'
##  is a treebase error that comes and goes....

## Libraries
library (treebase)
source (file.path ('tools', 'download_tools.R'))

## Parameters
if (!exists ('min.taxa')) {
  min.taxa <- 100
  subsample <- 100
  overwrite <- FALSE
}

## Dirs
output.dir <- file.path ('data', 'raw_trees', 'treebase')
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Setup
# First clear out all treebase trees from dir
treebase.files <- list.files (path = output.dir,
                              pattern = '\\.tre')
if (overwrite) {
  file.remove (file.path (output.dir, treebase.files))
  treebase.files <- NULL
}
# setup download log
download.log <- file.path (output.dir, 'metadata.csv')
if (file.exists (download.log)) {
  file.remove (download.log)
}
headers <- data.frame ("filename", "Study.id", "Tree.id", "kind",
                       "type", "quality", "ntaxa", "date", "publisher",
                       "author", "title" )
write.table (headers, download.log, sep = ',', row.names = FALSE,
             col.names = FALSE)

## Search
cat ('\nSearching suitable trees ....')
meta <- safeConnect (expr = suppressWarnings (metadata ()))
# choose those with more than min.taxa taxa
meta <- meta[meta$ntaxa >= min.taxa, ]
# choose those that are species trees
meta <- meta[meta$kind == 'Species Tree', ]
cat (paste0 ('\n[', nrow (meta), '] matching trees in TreeBase ....'))
cat (paste0 ('\n.... min size [', min (meta$ntaxa), ']'))
cat (paste0 ('\n.... mean size [', round (mean (meta$ntaxa)), ']'))
cat (paste0 ('\n.... max size [', max (meta$ntaxa), ']'))
# if subsample, choose a random set of trees in meta
if (is.numeric (subsample)) {
  if (subsample < nrow (meta)) {
    cat (paste0 ('\n.... taking [', subsample, '] subsample'))
    meta <- meta[sample (1:nrow (meta), subsample), ]
  }
}

## Download
cat ('\nDownloading trees ....')
counter <- 0
for (i in 1:nrow (meta)) {
  cat (paste0 ('\n.... tree [',i,']\n'))
  tree.data <- meta[i, ]
  filename <- paste0 (tree.data$Study.id, '.tre')
  # only download tree if it hasn't already been downloaded
  if (filename %in% treebase.files) {
    next
  }
  tree <- safeConnect (expr = {
    suppressWarnings (search_treebase (tree.data$Tree.id,
                                       by = 'id.tree',
                                       verbose = FALSE))})
  closeAllConnections ()
  if (length (tree) > 0) {
    tree.data <- cbind (filename, tree.data)
    # control for malformed downloaded trees
    error <- try (write.tree (
      tree[[1]], file = file.path (output.dir, filename),
      append = TRUE), silent = TRUE)
    if (!is.null (error)) {
      next
    }
    write.table (tree.data, download.log, sep = ',',
                 append = TRUE, col.names = FALSE,
                 row.names = FALSE)
    counter <- counter + 1
  } else {
    cat ('\n........ no tree could be retrieved')
  }
}
cat (paste0 ('\nStage complete, downloaded [', counter,'] trees'))