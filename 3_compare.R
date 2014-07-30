## 16/07/2014
## D.J. Bennett
## Comparing different model runs and natural trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Dirs
metadata <- read.csv (runlog, stringsAsFactors = FALSE)
data.dir <- 'data'

## Input
# Natural Trees
natural.trees <- list ()
tree.files <- list.files (data.dir, '.tre')
for (i in 1:length (tree.files)) {
  tree <- read.tree (file.path (data.dir, tree.files[i]))
  # choose the first tree if list
  if (class (tree) == 'multiPhylo') {
    tree <- tree[[1]]
  }
  natural.trees <- c (natural.trees, list (tree))
}
# Simulated trees
simulated.trees <- list ()
for (i in 1:nrow (metadata)) {
  res.dir <- file.path ('results', metadata$res.dir[i])
  tree <- read.tree (file.path (res.dir, 'MRMM.tre'))
  # use last tree of the model run
  tree <- tree[[length (tree)]]
  simulated.trees <- c (simulated.trees, list (tree))
}

## Calculate tree stats
# natural trees first
natural.tree.stats <- calcTreeShapeStats (natural.trees)

# calculate for each unique strength
simulated.tree.stats <- list ()
strengths <- sort (unique (metadata$strength))
for (i in 1:length (strengths)) {
  tree.i <- which (metadata$strength == strengths[i])
  temp.res <- list (calcTreeShapeStats (simulated.trees[tree.i]))
  names (temp.res) <- strengths[i]
  simulated.tree.stats <- c (simulated.tree.stats, temp.res)
}

## Plotting results
# extract the distribution of each stat and plot against ED strength
pdf (file.path ('results', 'treestats_ED_strength.pdf'))
stat.names <- c ('colless.stat', 'sackin.stat', 'iprime.stat',
                 'gamma.stat', 'tc.stat')
for (each in stat.names) {
  stats <- extractStat (simulated.tree.stats, each)
  x <- y <- NULL
  # reformat for a vector of x and y
  for (i in 1:length (stats)) {
    x <- c (x, rep (strengths[i], length (stats[[i]])))
    y <- c (y, stats[[i]])
  }
  plot (x = x, y = y, xlab = 'ED strength',
        ylab = paste0 ('Tree stat: [', each, ']'), pch = 19,
        col = rainbow (3, alpha = 0.8)[3])
  abline (lm (y ~ x))
  mtext (text = paste0 ('Mean natural stat: [',
                        mean (natural.tree.stats[[each]],
                              na.rm = TRUE), ']') )
}
closeDevices ()
# hist clade ages
#histClageAges (metadata)