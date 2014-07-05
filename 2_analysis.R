## 05/07/2014
## D.J. Bennett
## Analysis of modelled rise and fall of clades

## Parameters
res.dir <- 'ts200_int10_b1_d1_biasFP_date20140705'

## Dirs
res.dir <- file.path ('results', res.dir)
if (!file.exists (res.dir)) {
  stop (paste0 ('[', res.dir, '] does not exist!'))
}

## Libraries
source (file.path ('tools', 'analysis_tools.R'))

## Input
res <- read.csv (file = file.path (res.dir, 'clades_through_time.csv'))
time.steps <- res[ ,1]
res <- res[ ,-1]
trees <- read.tree (file = file.path (res.dir, 'ERMM.tre'))

## Generate figures
# plot clade successes across time
pdf (file.path (res.dir, 'clade_success.pdf'))
plotSuccess (res)
dev.off ()
# plot normalised clade success
pdf (file.path (res.dir, 'normalised_clade_success.pdf'))
plotNormalisedSuccess (res, min.time.span = 5, min.size = 1)
dev.off ()
# Create .gif of trees produced
plotTreeGrowth (trees, file.path (res.dir, 'ERMM_tree.gif'),
                time.steps)