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
res.dir <- paste0 ("t", time, 's', sample,
                   'b', round (birth), 'd', round (death),
                   'sn', seed.n, bias, '_',
                   format(Sys.time(), "%H%M_%a%b%d%Y"))
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
# calculate the adding and dropping of tips in order to
# know the number of iterations ahead of running model
add.bool.list <- calcAddBool (tree, birth, death, max.age = time,
                             sample.unit = sample)
# collect output for each interval
clade.performance <- list ()

## Model
runmodel <- function (i) {
  add.bool <- add.bool.list[[i]]
  # grow tree using ERMM
  tree <- growTree (add.bool, bias)
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
m_ply (.data = (i = 1:length (add.bool.list)), .fun = runmodel,
       .progress = 'time')
# convert list of dataframes into single dataframe
res <- reformat (clade.performance, sample)

## Saving results
write.csv (x = res, file = file.path (res.dir, 'clades_through_time.csv'))
