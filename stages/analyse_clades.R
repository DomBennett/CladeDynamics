## 05/07/2014
## D.J. Bennett
## Analysis of modelled MRMM trees
## This script is run through run.R

## TODO
# get this to work now with a multiPhylo
# import functions from model tools
# implement a new function for calculating CM and CG

# # calc 'success' of each node in tree
# temp.success <- .countChildren (tree, extinct)
# # add dataframe to a list
# clade.performance <<- c (clade.performance,
#                          list (temp.success))
# }
# # write seed tree first
# write.tree (tree, file = file.path (res.dir, 'MRMM.tre'))
# cat ('\nRunning tree growth model ...')
# m_ply (.data = (i = 1:iterations), .fun = runmodel,
#        .progress = 'time')
# # convert list of dataframes into single dataframe
# cat ('\nReformatting model output ...')
# res <- .reformat (clade.performance, sample)

## Dirs
# read in last run folder
res.dir <- read.csv (runlog, stringsAsFactors = FALSE)
res.dir <- res.dir[nrow (res.dir),'res.dir']
res.dir <- file.path ('results', as.character (res.dir))

## Libraries
source (file.path ('tools', 'analysis_tools.R'))

## Input
res <- read.csv (file = file.path (res.dir, 'clades_through_time.csv'))
time.steps <- res[ ,1]
res <- res[ ,-1]
trees <- read.tree (file = file.path (res.dir, 'MRMM.tre'))
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
closeDevices ()
# plot normalised clade success
cat ('\nPlotting normalised clade success ...')
pdf (file.path (res.dir, 'normalised_clade_success.pdf'))
plotNormalisedSuccess (res, min.time.span, min.size)
closeDevices ()
# Create .gif of trees produced
if (plot.tree.growth) {
  cat ('\nPlotting tree growth ...')
  plotTreeGrowth (trees, file.path (res.dir, 'MRMM_tree.gif'),
                  time.steps)
}
# plot fate ~ ED
cat ('\nPlotting Fate ~ ED ...')
pdf (file.path (res.dir, 'fate_ED.pdf'))
fates <- getFates (trees)
eds <- getEDs (trees)
plotFateVsED (fates, eds, time.lag = 1)
closeDevices ()