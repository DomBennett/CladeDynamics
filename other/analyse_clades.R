## 05/07/2014
## D.J. Bennett
## Analysis of modelled EDBMM tree clades

## Libraries
source (file.path ('tools', 'analysis_tools.R'))

## Input
trees <- read.tree (file.path (res.dir, treefilename))
# remove burnin, the first 10% of trees
first <- ceiling (length (trees)*0.1)
trees <- trees[first:length (trees)]

## Calculate
clades <- getCladeSuccess (trees, 1)
clade.stats <- calcCladeStats (clades)
# remove clades with time span < 5
clade.stats <- clade.stats[clade.stats$time.span > 5, ]
# remove all clades that do NOT appear and disappear
# i.e. I only want clades that come into being and go extinct
#  in this time span
clade.stats <- clade.stats[clade.stats$start > 1, ]
clade.stats <- clade.stats[clade.stats$end < nrow (clades), ]

## Save clade stats
filename <- sub ('\\.tre', '\\.csv', treefilename)
write.csv (x = clade.stats, file = file.path (res.dir, filename),
           row.names = FALSE)

## Generate figures
cat ('\nPlotting ...')
filename <- sub ('\\.tre', '\\.pdf', treefilename)
pdf (file.path (res.dir, filename))
plotSuccess (clades)
plotNormalisedSuccess (clades, 5, 5)
# Create .gif of trees produced
# if (plot.tree.growth) {
#   cat ('\nPlotting tree growth ...')
#   filename <- sub ('\\.tre', '\\.gif', treefilename)
#   plotTreeGrowth (trees, file.path (res.dir, filename),
#                   time.steps)
# }
# plot fate ~ ED (Won't work for trees without extinct)
# cat ('\nPlotting Fate ~ ED ...')
# pdf (file.path (res.dir, 'fate_ED.pdf'))
# fates <- getFates (trees)
# eds <- getEDs (trees)
# plotFateVsED (fates, eds, time.lag = 1)
closeDevices ()
