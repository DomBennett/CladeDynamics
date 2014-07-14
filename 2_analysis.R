## 05/07/2014
## D.J. Bennett
## Analysis of modelled rise and fall of clades

## Dirs
# read in last run folder
res.dir <- read.delim (file.path ('results', 'run_log.txt'),
                       header = FALSE, stringsAsFactors = FALSE)[ ,1]
#res.dir <- 'ts1000_int10_b1_d1_biasnone_date_seed2SunJul062014_time183238'
res.dir <- file.path ('results', res.dir [length (res.dir)])

## Libraries
source (file.path ('tools', 'analysis_tools.R'))

## Input
res <- read.csv (file = file.path (res.dir, 'clades_through_time.csv'))
time.steps <- res[ ,1]
res <- res[ ,-1]
trees <- read.tree (file = file.path (res.dir, 'ERMM.tre'))
# remove burnin
burnin.steps <- sum (burnin >= time.steps)
time.steps <- time.steps[-(1:burnin.steps)]
res <- res[ ,-(1:burnin.steps)]
trees <- trees[-(1:burnin.steps)]

## Generate figures
# plot clade successes across time
cat ('\nPlotting clade success ...')
pdf (file.path (res.dir, 'clade_success.pdf'))
plotSuccess (res)
dev.off ()
# plot normalised clade success
cat ('\nPlotting normalised clade success ...')
pdf (file.path (res.dir, 'normalised_clade_success.pdf'))
plotNormalisedSuccess (res, min.time.span, min.size)
dev.off ()
# Create .gif of trees produced
cat ('\nPlotting tree growth ...')
plotTreeGrowth (trees, file.path (res.dir, 'ERMM_tree.gif'),
                time.steps)
# plot fate ~ ED
cat ('\nPlotting Fate ~ ED ...')
pdf (file.path (res.dir, 'fate_ED.pdf'))
fates <- getFates (trees)
eds <- getEDs (trees)
plotFateVsED (fates, eds, time.lag = 1)
dev.off ()