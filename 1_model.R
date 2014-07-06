## 17/06/2014
## D.J. Bennett
## The rise and fall of clades: here I use
##  an equal rates markov model to add and remove
##  species to a tree. I then record the success
##  of each clade at set time intervals.
## N.B. It may crash unexpetedly as
##  all clades may have gone extinct; re-run

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Dirs
res.dir <- paste0 ("ts", time.steps, '_int', interval,
                   '_b', round (birth), '_d', round (death),
                   '_bias', bias, '_date', '_seed', seed.n,
                   format(Sys.time(), "%a%b%d%Y"),
                   '_time', format(Sys.time(), "%H%M%S"))
write.table (res.dir, file.path ('results', 'run_log.txt'),
             append = TRUE, row.names = FALSE, col.names = FALSE)
res.dir <- file.path ('results', res.dir)
if (!file.exists (res.dir)) {
  dir.create (res.dir)
}

## Globals
# count nodes and tips, necessary for adding new nodes
# otherwise would not be unique across time steps
max.node <- seed.n - 1
max.tip <- seed.n # starting tree has n tips and n-1 int node
extinct <- c () # vector of all extinct species

## Set-up
tree <- seedTree (seed.n)
# calc iterations, number of iterations each run of model
iterations <- time.steps/interval
# collect output for each interval
clade.performance <- list ()

## Model
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