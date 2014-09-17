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
stats[stats[['colless.stat']] > 1200,]

window.size <- 0.25
maxs <- seq (from = (-1.5 + window.size), to = 1,
             by = window.size)
mins <- seq (from = -1.5, to = (1 - window.size),
             by = window.size)
res <- rep (NA, length (maxs))
for (i in 1:length (maxs)) {
  temp <- stats[stats$strength < maxs[i] &
                  stats$strength > mins[i],
                'sackin.stat']
  indist <- rep (NA, 1000)
  for (j in 1:1000) {
    samp <- sample (real.stats, 1)
    indist[j] <- mean (temp) > samp
  }
  res[i] <- abs (0.5 - sum (indist)/1000)
}
par(xaxt="n")
plot (res, ylab = 'Difference from 0.5', pch = 19)
lablist <- paste0 (mins, ':', maxs)
axis(1, at = 1:length (maxs), labels = FALSE)
text(1:length (maxs), 0.01, labels = lablist, srt = 45, pos = 1, xpd = TRUE)

## Plotting results
# extract the distribution of each stat and plot against ED strength
cat ('\nPlotting ...')
pdf (file.path (res.dir, 'treestats_ED_strength.pdf'))
stat.names <- c ('colless.stat', 'sackin.stat', 'iprime.stat',
                 'gamma.stat', 'tc.stat')
for (each in stat.names) {
  x <- unlist (extractStat (tree.stats, each))
  y <- metadata$strength
  plot (x = x, y = y, ylab = 'Strength',
        xlab = paste0 ('Tree stat: [', each, ']'), pch = 19,
        col = rainbow (3, alpha = 0.8)[3])
  model <- lm (y ~ x)
  abline (model)
  #drawCorresPoints (model, stats[ ,each])
}
closeDevices ()