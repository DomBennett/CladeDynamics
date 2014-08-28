## 05/07/2014
## D.J. Bennett
## Analysis of modelled EDBMM tree clades

## TODO
# implement a new function for calculating CM and CG

## Libraries
source (file.path ('tools', 'analysis_tools.R'))

## Input
trees <- read.tree (file.path (res.dir, treefilename))
# remove burnin, the first 10% of trees
first <- ceiling (length (trees)*0.1)
trees <- trees[first:length (trees)]

## Calculate
clades <- getCladeSuccess (trees, 1)


## Generate figures
cat ('\nPlotting ...')
filename <- sub ('\\.tre', '\\.pdf', treefilename)
pdf (file.path (res.dir, filename))
plotSuccess (clades)
plotNormalisedSuccess (clades, 5, 5)
# Create .gif of trees produced
if (plot.tree.growth) {
  cat ('\nPlotting tree growth ...')
  filename <- sub ('\\.tre', '\\.gif', treefilename)
  plotTreeGrowth (trees, file.path (res.dir, filename),
                  time.steps)
}
# plot fate ~ ED (Won't work for trees without extinct)
# cat ('\nPlotting Fate ~ ED ...')
# pdf (file.path (res.dir, 'fate_ED.pdf'))
# fates <- getFates (trees)
# eds <- getEDs (trees)
# plotFateVsED (fates, eds, time.lag = 1)
closeDevices ()