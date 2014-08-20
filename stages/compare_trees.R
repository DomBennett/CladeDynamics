## 16/07/2014
## D.J. Bennett
## Comparing different model runs and natural trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Dirs
metadata <- read.csv (runlog, stringsAsFactors = FALSE)
data.dir <- 'data'

## Input
cat ('\nReading in data ...')
# load pre-calculated natural tree stats given leeway and target
filename <- paste0 ('natural_tree_stats_t', target, '_l', leeway, '.Rd')
if (!file.exists (file.path (data.dir, filename))) {
  stop ('No natural tree stats have been pre-calcualted with given parameters')
}
load (file.path (data.dir, filename))
# Read in last reconstructed phylogeny of simulated trees
trees <- list ()
for (i in 1:nrow (metadata)) {
  # read in tree
  tree.dir <- file.path (res.dir, metadata$treefilename[i])
  tree <- read.tree (tree.dir)
  # then add to list
  trees <- c (trees, list (tree))
}

## Calculate tree stats
cat ('\nCalculating tree stats ...')
tree.stats <- list ()
for (i in 1:length (trees)) {
  temp.res <- list (calcTreeShapeStats (trees[i]))
  tree.stats <- c (tree.stats, temp.res)
}

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
  drawCorresPoints (model, natural.tree.stats[ ,each])
}
closeDevices ()