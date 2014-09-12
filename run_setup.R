## 05/09/2014
## D.J. Bennett
## Download real trees from TreeBase and
##  calculate their tree shape stats
## Note. the results of this run should be
##  already be in /data

## Parameters
use.chronos <- FALSE # make trees ultrametric?
tree.dist <- 1 # the number of trees in a distribution for a polytomous tree
min.taxa <- 100 # the minimum tree size to be downloaded
targets <- c (100, 500, 1000) # the different tree sizes for which to calculate stats
leeway <- 10 # the percentage wobble around targets
iterations <- 100 # number of iterations of Yule comparison

## Stop if you don't want to run this!
x <- readline ('Running this script will delete already
           downloaded files, are you sure you want to continue?
           Hit return to continue. Ctrl+Z to exit.')
rm (x)

## Process
# download
cat ('\n--------------------------------')
cat (paste0 ('\n          Download'))
cat ('\n--------------------------------\n')
source (file.path ('stages', 'download.R'), print.eval = TRUE)
# parse
cat ('\n--------------------------------')
cat (paste0 ('\n          Parsing'))
cat ('\n--------------------------------\n')
source (file.path ('stages', 'parse.R'), print.eval = TRUE)
# precalculate
cat ('\n--------------------------------')
cat (paste0 ('\n          Precalculation'))
cat ('\n--------------------------------\n')
for (target in targets) {
  cat (paste0 ('\n--- Working on target [', target,'] ---'))
  source (file.path ('stages', 'precalculate.R'), print.eval = TRUE)
}
cat ('\n\nRun_setup.R complete')