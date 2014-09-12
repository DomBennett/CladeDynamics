## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for parsed trees

## Libraries
source (file.path ('tools', 'precalculate_tools.R'))
source (file.path ('tools', 'compare_tools.R'))

## Parameters
if (is.environment(.GlobalEnv)) {
  target <- 100
  leeway <- 10
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
min.n <- target - (target * leeway/100)
max.n <- target + (target * leeway/100)
treeinfo.master <- read.csv (file.path (input.dir,
                                        'treeinfo.csv'))
treeinfo <- data.frame ()
trees <- list ()
study.names <- NULL
for (i in 1:nrow (treeinfo.master)) {
  tree.file <- treeinfo.master[i,'filename']
  # read in
  tree <- read.tree (file.path (input.dir, tree.file))
  # pack into a multiphylo of right sized trees
  tree <- pack (tree, min.n = min.n, max.n = max.n)
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
for (set in trees) {
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
stats <- data.frame (colless.stat, sackin.stat, iprime.stat,
                                  gamma.stat, tc.stat)
stats <- cbind (treeinfo, stats)

## Save
filename <- paste0 ('t', target, '_l', leeway, '.Rd')
save (stats, file = file.path (output.dir, filename))
cat (paste0 ('\nStage complete, calculated for [', counter,'] tree sets'))