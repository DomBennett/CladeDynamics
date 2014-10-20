## 05/09/2014
## D.J. Bennett
## Download real trees from TreeBase and
##  calculate their tree shape stats
## Note. the results of this run should be
##  already be in /data

## Timestamp
cat (paste0 ('\nsetup.R started at [', Sys.time (), ']'))

## Parameters
use.chronos <- FALSE # make trees ultrametric?
subsample <- 10 # if numeric, take only subsample of trees for download
tree.dist <- 100 # the number of trees in a distribution for a polytomous tree
min.taxa <- 50 # the minimum tree size to be downloaded
max.taxa <- 200 # the maximum tree size to be downloaded
overwrite <- TRUE # delete all existing parsed trees and run again
reference <- TRUE # use a Yule reference for tree stats of same size

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
source (file.path ('stages', 'precalculate.R'), print.eval = TRUE)

## Timestamp
cat (paste0 ('\nsetup.R finished at [', Sys.time (), ']'))