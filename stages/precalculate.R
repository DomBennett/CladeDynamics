## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for parsed trees

## Libraries
source (file.path ('tools', 'precalculate_tools.R'))
source (file.path ('tools', 'compare_tools.R'))

## Parameters
if (!exists ('target')) {
  min.taxa <- 50
  max.taxa <- 200
  iterations <- 100
  reference  <- FALSE
}

## Dirs
input.dir <- file.path ('data', 'parsed_trees')
output.dir <- file.path ('data', 'treestats')
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Input
treeinfo.master <- read.csv (file.path (input.dir,
                                        'treeinfo.csv'))
treeinfo <- data.frame ()
trees <- list ()
study.names <- NULL
cat ('\nReading in trees and packing for target size ....')
for (i in 1:nrow (treeinfo.master)) {
  tree.file <- treeinfo.master[i,'filename']
  cat (paste0 ('\n.... working on [', tree.file,
               '] [', i, '/', nrow (treeinfo.master),']'))
  # read in
  tree <- read.tree (file.path (input.dir, tree.file))
  # pack into a multiphylo of right sized trees
  tree <- pack (tree, min.n = min.taxa, max.n = max.taxa)
  if (!is.null (tree)) {
    study.names <- c (study.names,
                      as.character (treeinfo.master[i,'Study.id']))
    # add to list
    trees <- c (trees, list (tree))
    # add its treeinfo
    treeinfo <- rbind (treeinfo, treeinfo.master[i, ])
  }
}
names (trees) <- study.names

## Calculate
counter <- 0
colless.stat <- sackin.stat <- iprime.stat <-
  gamma.stat <- tc.stat <- NULL
cat ('\nCalculating tree stats for sets of trees ....')
for (set in trees) {
  cat (paste0 ('\n.... working on set [', counter + 1,
               '/', length (trees),']'))
  stats <- calcTreeShapeStats (set, iterations = iterations,
                               reference = reference)
  # extract the mean value of the set
  colless.stat <- c (colless.stat, stats['mean.colless.stat'][[1]])
  sackin.stat <- c (sackin.stat, stats['mean.sackin.stat'][[1]])
  iprime.stat <- c (iprime.stat, stats['mean.iprime.stat'][[1]])
  gamma.stat <- c (gamma.stat, stats['mean.gamma.stat'][[1]])
  tc.stat <- c (tc.stat, stats['mean.tc.stat'][[1]])
  counter <- counter + 1
}
real.stats <- data.frame (colless.stat, sackin.stat, iprime.stat,
                                  gamma.stat, tc.stat)
real.stats <- cbind (treeinfo, real.stats)

## Save
filename <- paste0 ('min', min.taxa, '_max', max.taxa, '.Rd')
save (real.stats, file = file.path (output.dir, filename))
cat (paste0 ('\nStage complete, calculated for [', counter,'] tree sets'))