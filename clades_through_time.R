## 17/06/2014
## D.J. Bennett
## The rise and fall of clades: here I use
##  an equal rates markov model to add species to
##  a tree. I then record the success of each
##  clade at set time intervals.
## N.B. It may crash unexpetedly as
##  all clades may have gone extinct -- re-run
## TODO: record tree growth in a video/gif

## Libraries
library (ape)
library (MoreTreeTools)
library (plyr)
source ('functions.R')

## Parameters
time.steps <- 10000 # how many steps?
interval <- 10 # how often to record?
birth <- 3 # how many births to deaths?
death <- 1

## Globals
# count nodes, necessary for adding new nodes
# otherwise nodes would not be unique across time steps
max.node <- 1

## Set-up
# seed tree of two species
tree <- stree (2)
# add lengths and labels
tree$edge.length <- c (1,1)
tree$node.label <- 'n1'
# calc iterations, number of iterations each run of model
iterations <- time.steps/interval
# collect output for each interval
clade.performance <- list ()

## Model
runmodel <- function (iteration) {
  # grow tree using ERMM
  tree <- growTree (interval, birth, death)
  # calc 'success' of each node in tree
  current.success <- countChildren (tree)
  # add dataframe to a list
  clade.performance <<- c (clade.performance,
                          list (current.success))
}
m_ply (.data = (iteration = 1:iterations), .fun = runmodel,
       .progress = 'time')

## Recording results
# convert list of dataframes into single dataframe
res <- reformat (clade.performance, interval)
write.csv (x = res, file = 'clades_through_time.csv')
# plot clade successes across time
pdf ('clade_success.pdf')
plotSuccess (res)
dev.off ()
pdf ('normalised_clade_success.pdf')
plotNormalisedSuccess (res)
dev.off ()