## 16/07/2014
## D.J. Bennett
## Comparing different model runs and natural trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))

## Parameters
if (!exists ('res.dir')) {
  reference <- FALSE
  iterations <- 100
  target <- 100
  leeway <- 10
  res.dir <- file.path ('results', 'parameter_set_1')
  runlog <- file.path (res.dir, 'runlog.csv')
}

## Dirs
metadata <- read.csv (runlog, stringsAsFactors = FALSE)
data.dir <- file.path ('data', 'treestats')

## Input
cat ('\nReading in data ...')
# load pre-calculated natural tree stats given leeway and target
filename <- paste0 ('t', target, '_l', leeway, '.Rd')
if (!file.exists (file.path (data.dir, filename))) {
  stop ('No natural tree stats have been pre-calculated with given parameters')
}
load (file.path (data.dir, filename))
# Read in last reconstructed phylogeny of simulated trees
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

## Calculate tree stats
cat ('\nCalculating tree stats ...')
tree.stats <- list ()
for (i in 1:length (trees)) {
  temp.res <- list (calcTreeShapeStats (
    trees[i], iterations = iterations, reference = reference))
  tree.stats <- c (tree.stats, temp.res)
}

## Write out
save (tree.stats, file = file.path (
  res.dir, 'shapestats.Rd'))

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

res <- data.frame ()
window.size <- 0.5
stat.names <- c ('colless.stat', 'sackin.stat', 'iprime.stat',
                 'gamma.stat', 'tc.stat')
for (i in 1:length (stat.names)) {
  diff <- windowAnalysis (stats[ ,stat.names[i]],
                          real.stats[ ,stat.names[i]],
                          size = window.size,
                          strengths = stats[, 'strength'])
  temp <- data.frame (diff, size = window.size,
                      stat = stat.names[i])
  res <- rbind (res, temp)
}
res <- na.omit (res)

## Plotting results
cat ('\nPlotting ...')
pdf (file.path (res.dir, 'treestats_ED_strength.pdf'))
# plot distributions differences
p <- ggplot (res, aes (mid, diff))
p + geom_point (aes (colour = stat))
# plot stats against strength
for (each in stat.names) {
  x <- unlist (extractStat (tree.stats, each))
  y <- metadata$strength
  plot (x = x, y = y, ylab = 'Strength',
        xlab = paste0 ('Tree stat: [', each, ']'), pch = 19,
        col = rainbow (3, alpha = 0.8)[3])
  model <- lm (y ~ x)
  abline (model)
}
closeDevices ()