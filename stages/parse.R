## 11/09/2014
## D.J. Bennett
## Parse raw trees (lit and tb trees): make ultrametric and dichotomous

## Libraries
source (file.path ('tools', 'parse_tools.R'))

## Parameters
if (is.environment(.GlobalEnv)) {
  tree.dist <- 1 # how many trees in a dichotomous distribution?
  use.chronos <- FALSE # use chronos to make trees ultrametric?
}

## Dirs
treebase.dir <- file.path ('data', 'raw_trees', 'treebase')
literature.dir <- file.path ('data', 'raw_trees', 'literature')
output.dir <- file.path ('data', 'parsed_trees')
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}
# clear out all trees from dir
tree.files <- list.files (path = output.dir,
                              pattern = '\\.tre')
file.remove (file.path (output.dir, tree.files))

## Metadata
# read all tree data from lit and tb
metadata1 <- read.csv (file.path (treebase.dir, 'metadata.csv'))
tree.files <- file.path (treebase.dir, metadata1$filename)
metadata2 <- read.csv (file.path (literature.dir, 'metadata.csv'))
tree.files <- c (tree.files, file.path (literature.dir, metadata2$filename))
metadata <- rbind (metadata1, metadata2)
# choose the first study.id metadata if duplicated
#  (not the best solution but good enough)
tree.files <- tree.files[!duplicated (metadata['Study.id'])]
metadata <- metadata[!duplicated (metadata['Study.id']), ]
rm (metadata1, metadata2)

## Data slots
# add extra data slots to metadata, to track for 
#  causes of differences in tree shape
metadata$multi <- FALSE # part of a collection of trees
metadata$poly <- FALSE # original tree was polytomous
metadata$bl <- TRUE # original tree had branch lengths
metadata$ultra <- FALSE # original tree was ultrametric
metadata$chronos <- FALSE # made ultrametric by chronos

## Parse
counter <- 0
# container for multiPhylos
trees <- list ()
# metadata for trees that get used
treeinfo <- data.frame ()
for (i in 1:length (tree.files)) {
  tempinfo <- metadata[i, ]
  # suppress add.terminal warnings
  tree <- suppressWarnings (read.tree (tree.files[i]))
  # choose the biggest tree of the trees for a study
  if (class (tree) == 'multiPhylo') {
    sizes <- unlist(lapply (tree, getSize))
    tree <- tree[[which (sizes == max(sizes))[1]]]
    tempinfo['multi'] <- TRUE
  }
  # does it have branch lengths?
  bl.bool <- !is.null (tree$edge.length) &&
    all (!is.na (tree$edge.length))
  if (!bl.bool) {
    # make sure if any edgelenths are NA,
    #  edgelenghts are removed
    tree$edge.length <- NULL
    tempinfo['bl'] <- FALSE
  }
  # is it ultrametric?
  ultra.bool <- bl.bool && is.ultrametric (tree)
  if (ultra.bool) {
    tempinfo['ultra'] <- TRUE
  }
  # is it polytomous?
  poly.bool <- getSize (tree) != (tree$Nnode + 1)
  # print progress
  cat (paste0 ('\nWorking on [', tempinfo['Tree.id'],']'))
  # if not ultrametric make it (if I can)
  if (use.chronos && bl.bool && !ultra.bool) {
    cat ('\n.... using chronos')
    tree <- safeChronos (tree, lambda = 1, quiet = TRUE)
    if (is.ultrametric (tree)) {
      class (tree) <- 'phylo'
      tempinfo['chronos'] <- TRUE
    }
  }
  # if polytomous, convert to a distribution
  if (poly.bool) {
    tempinfo['poly'] <- TRUE
    tree <- convertToDist (tree)
  }
  # convert to multiPhylo before adding to trees
  if (any (class (tree) != 'multiPhylo')) {
    tree <- list (tree)
    class (tree) <- 'multiPhylo'
  }
  trees <- c (trees, list (tree))
  treeinfo <- rbind (treeinfo, tempinfo)
  counter <- counter + 1
}

## Output
for (i in 1:length (trees)) {
  write.tree (trees[[i]], file.path (
    output.dir, treeinfo[i, 'filename']))
}
write.csv (treeinfo, file.path (output.dir, 'treeinfo.csv'),
           row.names = FALSE)
cat (paste0 ('\nStage complete, parsed [', counter,'] trees'))