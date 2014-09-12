## 07/08/2014
## D.J. Bennett
## Model trees using EDBMM

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Parameters
if (is.environment(.GlobalEnv)) {
  # run outside of run.R for testing stage
  res.dir <- file.path ('results', 'test_parameters')
  list.files (res.dir)
  runlog <- file.path (res.dir, 'runlog.csv')
  if (file.exists (runlog)) {
    file.remove (runlog)
  }
  headers <- data.frame ('treefilename', 'strength', 'bias',
                         'birth', 'death')
  write.table (headers, runlog, sep = ',', row.names = FALSE,
               col.names = FALSE)
  min.n <- 40
  max.n <- 60
  min.strength <- -1
  max.strength <- 1
  birth <- 2
  death <- 1
  bias <- 'FP'
  stop.by <- 'n'
  seed <- 2
  n <- 5
}


## Loop
for (j in 1:n) {
  # Setup
  cat (paste0 ('\n--- Working on model [', j,'] of [', n,
               '] ---'))
  strength <- runif (1, min.strength, max.strength)
  stop.at <- round (runif (1, min.n, max.n))
  
  # Dirs
  counter <- nrow (read.csv (runlog))
  # create a tree filename
  treefilename <- paste0 ('tree', counter, '.tre')
  # record parameters in runlog
  parameters <- data.frame (treefilename, strength, bias,
                            birth, death)
  write.table (parameters, runlog, sep = ',', append = TRUE,
               col.names = FALSE, row.names = FALSE)
  
  # Model
  cat (paste0 ('\nModelling tree of size [', stop.at, '(',
               stop.by, ')] ...'))
  # if seed is greater than 2
  if (seed > 2) {
    # ... create a seed tree of 2-1
    seed.tree <- runEDBMM (birth = 2, death = 1,
                           stop.at = seed, stop.by = 'n',
                           strength = strength, bias = bias,
                           fossils = FALSE, record = FALSE)
    # reset names so no overlap
    seed.tree$tip.label <- paste0 ('t', 1:getSize (seed.tree))
    seed.tree$node.label <- paste0 ('n', 1:seed.tree$Nnode)
    tree <- runEDBMM (birth = birth, death = death,
                      stop.at = stop.at, stop.by = stop.by,
                      strength = strength, bias = bias,
                      fossils = FALSE, record = FALSE,
                      seed.tree = seed.tree)
  } else {
    tree <- runEDBMM (birth = birth, death = death,
                      stop.at = stop.at, stop.by = stop.by,
                      strength = strength, bias = bias,
                      fossils = FALSE, record = FALSE)
  }
  
  # Write out
  cat ('\nSaving ...')
  write.tree (tree, file = file.path (res.dir, treefilename))
}