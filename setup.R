## 05/09/2014
## D.J. Bennett
## Download real trees from TreeBase and
##  calculate their tree shape stats
## Note. the results of this run should be
##  already be in /data

## Timestamp
cat (paste0 ('\nsetup.R started at [', Sys.time (), ']'))

## Parameters
use.chronos <- TRUE # make trees ultrametric?
subsample <- 100 # if numeric, take only subsample of trees for download
tree.dist <- 100 # the number of trees in a distribution for a polytomous tree
min.taxa <- 100 # the minimum tree size to be downloaded
targets <- c (100, 200, 300) # the different tree sizes for which to calculate stats
leeway <- 10 # the percentage wobble around targets
overwrite <- FALSE # delete all existing parsed trees and run again

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

## Timestamp
cat (paste0 ('\nsetup.R finished at [', Sys.time (), ']'))