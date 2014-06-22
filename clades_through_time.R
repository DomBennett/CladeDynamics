## 17/06/2014
## D.J. Bennett
## The rise and fall of clades: here I use
##  an equal rates markov model to add and remove
##  species to a tree. I then record the success
##  of each clade at set time intervals.
## N.B. It may crash unexpetedly as
##  all clades may have gone extinct; re-run
## TODO: store tree with fossil species, for plotting
## TODO: I don't think the node allocation is working,
##  clades are disappearing far two quickly to be the
##  of chance.

## Libraries
library (ape)
library (MoreTreeTools)
library (plyr)
source ('functions.R')

## Dirs
figs.dir <- "figures"
if (!file.exists (figs.dir)) {
  dir.create (figs.dir)
}
res.dir <- "results"
if (!file.exists (res.dir)) {
  dir.create (res.dir)
}

## Parameters
time.steps <- 3000 # how many steps?
# ~63 mins for 10000 on MacBook Air 1.4Ghrz*2
interval <- 10 # how often to record?
birth <- 1.1 # how many births to deaths?
death <- 1
pe.bias <- TRUE # Bias by PE?

## Globals
# count nodes and tips, necessary for adding new nodes
# otherwise would not be unique across time steps
max.node <- 1
max.tip <- 2 # starting tree has 2 tips and 1 int node

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
  tree <- growTree (interval, birth, death, pe.bias)
  # write tree to disk
  write.tree (tree, file = file.path (
    res.dir, 'ERMM.tre'), append = TRUE)
  # calc 'success' of each node in tree
  current.success <- countChildren (tree)
  # add dataframe to a list
  clade.performance <<- c (clade.performance,
                          list (current.success))
}
if (file.exists (file.path (res.dir, 'ERMM.tre'))) {
  file.remove (file.path (res.dir, 'ERMM.tre'))
}
m_ply (.data = (iteration = 1:iterations), .fun = runmodel,
       .progress = 'time')

## Saving results
# convert list of dataframes into single dataframe
res <- reformat (clade.performance, interval)
write.csv (x = res, file = file.path (res.dir, 'clades_through_time.csv'))
#res <- read.csv (file = file.path (res.dir, 'clades_through_time.csv'))
# plot clade successes across time
pdf (file.path (figs.dir, 'clade_success.pdf'))
plotSuccess (res)
dev.off ()
pdf (file.path (figs.dir, 'normalised_clade_success.pdf'))
plotNormalisedSuccess (res, min.time.span = 30)
dev.off ()
# Create .gif of trees produced
# read trees
trees <- read.tree (file = file.path (res.dir, 'ERMM.tre'))
time.step <- interval
par (mar = c (2, 2, 3, 2) + 0.1)
png (file = "temp_tree_%010d.png")
for (each in trees) {
  plot (each, show.tip.label = FALSE)
  mtext (paste0 ('t = ', time.step), adj = 0.1, line = 0)
  time.step <- time.step + interval
}
dev.off()
system (paste0 ("convert -delay 10 *.png ",
                file.path (figs.dir, 'ERMM_tree.gif')))
file.remove (list.files (pattern = ".png"))