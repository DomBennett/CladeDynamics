calcTreeStats <- function (trees) {
  engine <- function (i) {
    tree <- trees[[i]]
    # imbalance stats
    #ts.tree <- as.treeshape (tree)
    #colless.stat <- colless (ts.tree, 'yule')
    #sackin.stat <- sackin (ts.tree, 'yule')
    # branching stats
    if (is.null (tree$edge.length)) {
      gamma.stat <- psv.stat <- age <- pd <- NA
    } else {
      # get gamma
      gamma.stat <- gammaStat (tree)
      # to get PSV, create community matrix
      samp <- matrix (rep (1, getSize (tree) * 2), nrow=2)
      colnames (samp) <- tree$tip.label
      psv.res <- psv (samp, tree)
      psv.stat <- psv.res[1,1]
      # get age
      age <- getSize (tree, 'rtt')
      # get pd
      pd <- getSize (tree, 'pd')
    }
    data.frame (gamma = gamma.stat, psv = psv.stat, age, pd)
  }
  res_part <- mdply (.data = data.frame (i = 1:length (trees)), .fun = engine)[, -1]
}


nplys <- sapply(trees, function(t) sum(t$edge.length == 0))
res_part$nply <- nplys
plot(psv~nply, data=res_part)

## Libraries
source (file.path ('tools', 'precalculate_tools.R'))

if (!exists ('min.taxa')) {
  min.taxa <- 50
  max.taxa <- 500
  rate.smooth <- 'pathD8'
}

## Dirs
input.dir <- file.path ('data', paste0 ('parsed_trees_', rate.smooth))
output.dir <- file.path ('data', 'treestats')
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

treeinfo.master <- read.csv (file.path (input.dir,
                                        'treeinfo.csv'))
treeinfo.master <- treeinfo.master[treeinfo.master$quality == 'literature', ]
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
  tree <- read.tree (tree.file)
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
colless <- sackin <- psv <- gamma <- age <- pd <- rep (NA, length (trees))
cat ('\nCalculating tree stats for sets of trees ....')
for (i in 1:length (trees)) {
  set <- trees[[i]]
  cat (paste0 ('\n.... working on set [', i,
               '/', length (trees),']'))
  stats <- try (expr= {calcTreeStats (set)}, silent = TRUE)
  if (class (stats) == 'try-error') {
    cat (paste0 ('\n.... skipping [', counter+1, '] the following error was encountered:\n',
                 attr(stats, 'condition')))
    next
  }
  # extract mean stats of the set
  colless[i] <- mean (stats[ ,'colless'], na.rm = TRUE)
  sackin[i] <- mean (stats[ ,'sackin'], na.rm = TRUE)
  psv[i] <- mean (stats[ ,'psv'], na.rm = TRUE)
  gamma[i] <- mean (stats[ ,'gamma'], na.rm = TRUE)
  age[i] <- mean (stats[ ,'age'], na.rm = TRUE)
  pd[i] <- mean (stats[ ,'pd'], na.rm = TRUE)
  counter <- counter + 1
}
real.stats <- data.frame (colless, sackin, psv, gamma, age, pd)
real.stats <- cbind (treeinfo, real.stats)

lngths <- unlist(lapply(tree, function(t) length(t$tip.label)))
sum(lngths)



# ADD MISSING DATA
source (file.path ('tools', 'compare_tools.R'))
data.dir <- file.path ('data', 'treestats')
empirical.filenames <- c ('min50_max500.Rd',
                          'min50_max500_rspathD8.Rd',
                          'min50_max500_rschronopl.Rd',
                          'min50_max500_rschronoMPL.Rd')
real.stats <- readInMultiple (file.path (data.dir, empirical.filenames))
# record gravity metrics for originally ultrametric trees
ultra_psv <- real.stats$psv.[real.stats$ultra]
ultra_gamma <- real.stats$gamma.[real.stats$ultra]
real.stats$gamma. <- real.stats$psv. <- NULL

# add ultrametrics to each of the rs results
real.stats$psv.pathD8[real.stats$ultra] <- ultra_psv
real.stats$psv.chronopl[real.stats$ultra] <- ultra_psv
real.stats$psv.chronoMPL[real.stats$ultra] <- ultra_psv
real.stats$gamma.pathD8[real.stats$ultra] <- ultra_gamma
real.stats$gamma.chronopl[real.stats$ultra] <- ultra_gamma
real.stats$gamma.chronoMPL[real.stats$ultra] <- ultra_gamma