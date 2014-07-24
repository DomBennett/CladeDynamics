## 17/06/2014
## D.J. Bennett
## Modelling tree growth using a modified rates
## markov model, where rates are determined by ED.
## This script is run through run.R

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Dirs
# create a unique dir based on time for every run
counter <- nrow (read.csv (runlog))
res.dir <- paste0 ('n', counter, '_',
                   format (Sys.time (), "%H%M_%d%m%y"))
# record parameters in runlog
parameters <- data.frame (res.dir, strength, bias, time,
                          sample, birth, death, seed.n,
                          min.time.span, min.size)
write.table (parameters, runlog, sep = ',', append = TRUE,
             col.names = FALSE, row.names = FALSE)
# set res.dir's path
res.dir <- file.path ('results', res.dir)
if (!file.exists (res.dir)) {
  dir.create (res.dir)
}

## Set-up
# globals: count nodes and tips, necessary for adding new nodes
#  otherwise would not be unique across time steps
max.node <- seed.n - 1
max.tip <- seed.n # starting tree has n tips and n-1 int node
extinct <- c () # vector of all extinct species
# seed tree age is the same as a sample unit
tree <- seedTree (seed.n, sample)
# collect output for each sample point
clade.performance <- list ()
# calculate iterations
iterations <- time/sample

## Model
runmodel <- function (i) {
  # grow tree using MRMM
  tree <<- growMRMMTree (birth = birth, death = death,
                        stop.at = sample, stop.by = 'max.time',
                        strength = strength, bias = bias,
                        seed.tree = tree)
  # write last tree to disk
  write.tree (tree, file = file.path (
    res.dir, 'MRMM.tre'), append = TRUE)
  # calc 'success' of each node in tree
  temp.success <- .countChildren (tree, extinct)
  # add dataframe to a list
  clade.performance <<- c (clade.performance,
                          list (temp.success))
}
# write seed tree first
write.tree (tree, file = file.path (res.dir, 'MRMM.tre'))
cat ('\nRunning tree growth model ...')
m_ply (.data = (i = 1:iterations), .fun = runmodel,
       .progress = 'time')
# convert list of dataframes into single dataframe
cat ('\nReformatting model output ...')
res <- .reformat (clade.performance, sample)

## Saving results
cat ('\nSaving results ...')
write.csv (x = res, file = file.path (res.dir, 'clades_through_time.csv'))