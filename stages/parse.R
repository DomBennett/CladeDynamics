## 11/09/2014
## D.J. Bennett
## Parse raw trees (lit and tb trees): make ultrametric and dichotomous

## Libraries
library (foreach)
library (doMC)
source (file.path ('tools', 'parse_tools.R'))

## Functions
getTreeFiles <- function (dirs) {
  ## Add new tree files and metadata to metadata and
  ##  tree.files from vector of dirs
  for (d in dirs) {
    if (file.exists (file.path (d, 'metadata.csv'))) {
      tmetadata <- read.csv (file.path (d, 'metadata.csv'))
      bool <- !tmetadata$filename %in% deja.vues
      tree.files <<- c (tree.files, file.path (
        d, tmetadata$filename[bool]))
      metadata <<- rbind (metadata, tmetadata[bool, ])
    }
  }
}

## Parameters
if (!exists ('tree.dist')) {
  tree.dist <- 1 # how many trees in a dichotomous distribution?
  rate.smooth <- 'chronopl' # none, pathD8, chronoMPL or chronopl
  overwrite <- TRUE
  subsample <- FALSE
  ncpus <- 2
}
registerDoMC (ncpus)

## Dirs
treebase.dir <- file.path ('data', 'raw_trees',
                           'treebase')
literature.dir <- file.path ('data', 'raw_trees',
                             'literature')
output.dir <- file.path ('data', paste0 ('parsed_trees_', rate.smooth))
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}
parse.log <- file.path (output.dir, 'treeinfo.csv')

## Overwrite
if (overwrite) {
  # clear out all trees from dir
  tree.files <- list.files (path = output.dir,
                            pattern = '\\.tre')
  file.remove (file.path (output.dir, tree.files))
  if (file.exists (parse.log)) {
    file.remove (parse.log)
  }
  # metadata for trees that get used
  treeinfo <- data.frame ()
  deja.vues <- NULL
} else {
  deja.vues <- list.files(output.dir, pattern = '\\.tre')
}
if (!file.exists (parse.log)) {
  headers <- data.frame ("filename", "Study.id",
                         "Tree.id", "kind",
                         "type", "quality",
                         "ntaxa", "date",
                         "publisher", "author",
                         "title", "multi",
                         "poly", "bl", "ultra",
                         "rate.smooth")
  write.table (headers, parse.log, sep = ',',
               row.names = FALSE, col.names = FALSE)
}

## Metadata + files
# read all tree data from lit and tb
tree.files <- metadata <- NULL
getTreeFiles (c (treebase.dir, literature.dir))
if (length (tree.files) == 0) {
  stop ('Stopped: Overwrite false and no new files to add!')
}
# choose the first study.id metadata if duplicated
#  (not the best solution but good enough)
tree.files <- tree.files[!duplicated (metadata['Study.id'])]
metadata <- metadata[!duplicated (metadata['Study.id']), ]
# if more tree.files than subsample, choose subsample at random
if (is.numeric (subsample)) {
  if (subsample < length (tree.files)) {
    cat (paste0 ('.... taking [', subsample, '] subsample'))
    random.is <- sample (1:length (tree.files), subsample)
    tree.files <- tree.files[random.is]
    metadata <- metadata[random.is, ]
  }
}

## Data slots
# add extra data slots to metadata, to track for 
#  causes of differences in tree shape
metadata$multi <- FALSE # part of a collection of trees
metadata$poly <- FALSE # original tree was polytomous
metadata$bl <- TRUE # original tree had branch lengths
metadata$ultra <- FALSE # original tree was ultrametric
metadata$rate.smooth <- FALSE # made ultrametric by rate.smooth

## Parse
counter <- foreach (i=1:length (tree.files)) %dopar% {
  tempinfo <- metadata[i, ]
  # suppress add.terminal warnings
  tree <- suppressWarnings (read.tree (tree.files[i]))
  # choose the biggest tree of the trees for a study
  if (class (tree) == 'multiPhylo') {
    sizes <- unlist(lapply (tree, getSize))
    tree <- tree[[which (sizes == max(sizes))[1]]]
    tempinfo['multi'] <- TRUE
  }
  # ensure we have a tree
  if (class (tree) == 'phylo') {
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
    cat (paste0 ('\nWorking on [', tree.files[i],
                 '] [', i, '/', length (tree.files),']'))
    # if not ultrametric make it (if I can)
    if (rate.smooth != 'none' && bl.bool && !ultra.bool) {
      rs.tree <- runRateSmoother (tree, i=i, rsmoother=rate.smooth)
      if ('phylo' %in% class (rs.tree) && is.ultrametric (rs.tree)) {
        class (rs.tree) <- 'phylo'
        tree <- rs.tree
        tempinfo['rate.smooth'] <- rate.smooth
      }
    }
    # if polytomous, convert to a distribution
    if (poly.bool) {
      cat ('\n.... converting to distribution')
      tempinfo['poly'] <- TRUE
      tree <- try (convertToDist (tree), silent = TRUE)
    }
    if (class (tree) == 'try-error') {
      cat (paste0 ('\n.... error occurred: [', tree, ']'))
      cat ('....\n moving to next tree')
    } else {
      # convert to multiPhylo before adding to trees
      if (any (class (tree) != 'multiPhylo')) {
        tree <- list (tree)
        class (tree) <- 'multiPhylo'
      }
      # write out
      write.tree (tree, file.path (
        output.dir, tempinfo[['filename']]))
      # save details
      write.table (tempinfo, parse.log, sep = ',',
                   append = TRUE, col.names = FALSE,
                   row.names = FALSE)
      1
    }
  }
}
counter <- sum (unlist (counter))
cat (paste0 ('\nStage complete, parsed [', counter,'] trees'))