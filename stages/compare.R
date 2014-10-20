## 16/07/2014
## D.J. Bennett
## Comparing different model runs and natural trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))

## Parameters
if (!exists ('analysis.parameters')) {
  testset <- list (n.model = 2, seed = 2,
                   max.birth = 5, min.birth = 1.1,
                   max.death = 1, min.death = 1,
                   bias = 'FP', stop.by = 'n',
                   max.ntaxa = 200, min.ntaxa = 50,
                   min.psi = -1, max.psi = 1,
                   reference = TRUE, iterations = 100)
  analysis.parameters <- list (testset = testset)
  rm (testset)
}

## Dirs
data.dir <- file.path ('data', 'treestats')

## Functions
readTrees <- function (metadata, res.dir, runlog) {
  trees <- list ()
  for (i in 1:nrow (metadata)) {
    # read in tree
    tree.dir <- file.path (res.dir, metadata$treefilename[i])
    tree <- read.tree (tree.dir)
    # make sure it isn't multiPhylo
    if (class (tree) == 'multiPhylo') {
      tree <- tree[[length (tree)]]
    }
    # then add to list
    trees <- c (trees, list (tree))
  }
  trees
}

calcStats <- function (trees, metadata, iterations, reference,
                       res.dir) {
  # loop through each tree
  tree.stats <- list ()
  for (i in 1:length (trees)) {
    temp.res <- list (calcTreeShapeStats (
      trees[i], iterations = iterations, reference = reference))
    tree.stats <- c (tree.stats, temp.res)
  }
  # write out
  save (tree.stats, file = file.path (res.dir, 'shapestats.Rd'))
  # TODO: rethink how the calcTreeShapeStats works to avoid this loop
  colless.stat <- sackin.stat <- iprime.stat <- 
    gamma.stat <- tc.stat <- rep (NA, length (tree.stats))
  for (i in 1:length (tree.stats)) {
    colless.stat[i] <- tree.stats[[i]][['colless.stat']]
    sackin.stat[i] <- tree.stats[[i]][['sackin.stat']]
    iprime.stat[i] <- tree.stats[[i]][['iprime.stat']]
    gamma.stat[i] <- tree.stats[[i]][['gamma.stat']]
    tc.stat[i] <- tree.stats[[i]][['tc.stat']]
  }
  stats <- data.frame (colless.stat, sackin.stat, iprime.stat,
                       gamma.stat, tc.stat)
  stats <- cbind (metadata, stats)
  stats
}

compare <- function (stats, real.stats) {
  res <- data.frame ()
  window.size <- 0.5
  stat.names <- c ('colless.stat', 'sackin.stat', 'iprime.stat',
                   'gamma.stat', 'tc.stat')
  for (i in 1:length (stat.names)) {
    diff <- windowAnalysis (stats[ ,stat.names[i]],
                            real.stats[ ,stat.names[i]],
                            size = window.size,
                            psi = stats[, 'psi'], incr = 0.1)
    temp <- data.frame (diff, size = window.size,
                        stat = stat.names[i])
    res <- rbind (res, temp)
  }
  res <- na.omit (res)
  res
}

plotResults <- function (res, stats, res.dir) {
  pdf (file.path (res.dir, 'treestats_ED_strength.pdf'))
  # plot distributions differences
  p <- ggplot (res, aes (mid, diff))
  print (p + geom_point (aes (colour = stat)))
  # plot stats against strength
  stat.names <- c ('colless.stat', 'sackin.stat', 'iprime.stat',
                   'gamma.stat', 'tc.stat')
  for (each in stat.names) {
    x <- stats[ ,each]
    y <- metadata$psi
    plot (x = x, y = y, ylab = expression (psi),
          xlab = paste0 ('Tree stat: [', each, ']'), pch = 19,
          col = rainbow (3, alpha = 0.8)[3])
    model <- lm (y ~ x)
    abline (model)
  }
  closeDevices ()
}

iterateAnalysis <- function (i) {
  ## Iterate through analyses
  cat ('\n--------------------------------')
  cat (paste0 ('\n          Analysis [', i, ']'))
  cat ('\n--------------------------------\n')
  # get res.dir and runlog
  analysis.name <- names (analysis.parameters)[i]
  res.dir <- file.path ('results', analysis.name)
  runlog <- file.path (res.dir, 'runlog.csv')
  pars <- analysis.parameters[[i]]
  cat ('\nReading in data ...')
  # get metadata
  metadata <- read.csv (runlog, stringsAsFactors = FALSE)
  # load pre-calculated natural tree stats
  filename <- paste0 ('min', pars$min.ntaxa, '_max', pars$max.ntaxa, '.Rd')
  load (file.path (data.dir, filename))
  trees <- readTrees (metadata, res.dir, runlog)
  cat ('\nCalculating tree stats ...')
  stats <- calcStats (trees, metadata, pars$iterations, pars$reference, res.dir)
  cat ('\nCompare to real trees ...')
  res <- compare (stats, real.stats)
  cat ('\nPlotting ...')
  plotResults (res, stats, res.dir)
}

## Run
m_ply (.data = data.frame (i = 1:length (analysis.parameters)),
       .fun = iterateAnalysis)