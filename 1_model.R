## 17/06/2014
## D.J. Bennett
## The rise and fall of clades: here I use
##  an equal rates markov model to add and remove
##  species to a tree. I then record the success
##  of each clade at set time intervals.
## N.B. It may crash unexpetedly as
##  all clades may have gone extinct; re-run

## Parameters
time.steps <- 200 # how many steps?
interval <- 10 # how often to record?
birth <- 1.1 # how many births to deaths per unit of branch length?
death <- 1
bias <- 'FP' # 'none', 'PE' or 'FP'

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Dirs
res.dir <- paste0 ("ts", time.steps, '_int', interval,
                   '_b', round (birth), '_d', round (death),
                   '_bias', bias, '_date', gsub ('-', '', Sys.Date ()))
res.dir <- file.path ('results', res.dir)
if (!file.exists (res.dir)) {
  dir.create (res.dir)
}

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
extinct <- c () # vector of all extinct species
runmodel <- function (iteration) {
  # grow tree using ERMM
  tree <- growTree (interval, birth, death, bias)
  # write tree to disk
  write.tree (tree, file = file.path (
    res.dir, 'ERMM.tre'), append = TRUE)
  # calc 'success' of each node in tree
  temp.success <- countChildren (tree)
  # add dataframe to a list
  clade.performance <<- c (clade.performance,
                          list (temp.success))
}
if (file.exists (file.path (res.dir, 'ERMM.tre'))) {
  file.remove (file.path (res.dir, 'ERMM.tre'))
}
m_ply (.data = (iteration = 1:iterations), .fun = runmodel,
       .progress = 'time')
# convert list of dataframes into single dataframe
res <- reformat (clade.performance, interval)

## Saving results
write.csv (x = res, file = file.path (res.dir, 'clades_through_time.csv'))