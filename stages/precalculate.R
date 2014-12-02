## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for parsed trees

## Libraries
source (file.path ('tools', 'precalculate_tools.R'))
source (file.path ('tools', 'compare_tools.R'))

## Parameters
if (!exists ('min.taxa')) {
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
  tree.file <- file.path (input.dir, treeinfo.master[i,'filename'])
  if (!file.exists (tree.file)) {
    next
  }
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
colless <- sackin <- gamma <- tci <- NULL
cat ('\nCalculating tree stats for sets of trees ....')
for (set in trees) {
  cat (paste0 ('\n.... working on set [', counter + 1,
               '/', length (trees),']'))
  stats <- try (expr= {calcTreeStats (set)}, silent = TRUE)
  if (class (stats) == 'try-error') {
    cat (paste0 ('\n.... skipping [', counter+1, '] the following error was encountered:\n', attr(stats, 'condition')))
    next
  }
  # extract mean stats of the set
  colless <- c (colless, mean (stats[ ,'colless'], na.rm = TRUE))
  sackin <- c (sackin, mean (stats[ ,'sackin'], na.rm = TRUE))
  gamma <- c (gamma, mean (stats[ ,'gamma'], na.rm = TRUE))
  tci <- c (tci, mean (stats[ ,'tci'], na.rm = TRUE))
  counter <- counter + 1
}
real.stats <- data.frame (colless, sackin, gamma, tci)
real.stats <- cbind (treeinfo, real.stats)

## Save
filename <- paste0 ('min', min.taxa, '_max', max.taxa, '.Rd')
save (real.stats, file = file.path (output.dir, filename))
cat (paste0 ('\nStage complete, calculated for [', counter,'] tree sets'))