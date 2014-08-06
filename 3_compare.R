## 16/07/2014
## D.J. Bennett
## Comparing different model runs and natural trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Dirs
metadata <- read.csv (runlog, stringsAsFactors = FALSE)
data.dir <- 'data'

## Input
# load pre-calculated natural tree stats
load (file.path (data.dir, 'natural_tree_stats.Rd'))
# Read in last reconstructed phylogeny of simulated trees
simulated.trees <- list ()
for (i in 1:nrow (metadata)) {
  res.dir <- file.path ('results', metadata$res.dir[i])
  tree <- read.tree (file.path (res.dir, 'MRMM.tre'))
  # use last tree of the model run
  tree <- tree[[length (tree)]]
  # drop all extinct nodes
  tree <- drop.extinct (tree)
  # then add to list
  simulated.trees <- c (simulated.trees, list (tree))
}

## Calculate tree stats
simulated.tree.stats <- list ()
for (i in 1:length (simulated.trees)) {
  temp.res <- list (calcTreeShapeStats (simulated.trees[i]))
  simulated.tree.stats <- c (simulated.tree.stats, temp.res)
}

## Plotting results
# extract the distribution of each stat and plot against ED strength
pdf (file.path ('results', 'treestats_ED_strength.pdf'))
stat.names <- c ('colless.stat', 'sackin.stat', 'iprime.stat',
                 'gamma.stat', 'tc.stat')
for (each in stat.names) {
  y <- unlist (extractStat (simulated.tree.stats, each))
  x <- metadata$strength
  plot (x = x, y = y, xlab = 'ED strength',
        ylab = paste0 ('Tree stat: [', each, ']'), pch = 19,
        col = rainbow (3, alpha = 0.8)[3])
  model <- lm (y ~ x)
  Y <- mean (natural.tree.stats[[each]], na.rm = TRUE)
  X <- predict (model, newdata = data.frame (y = Y))
}
closeDevices ()
# hist clade ages
#histClageAges (metadata)