## 07/08/2014
## D.J. Bennett
## Run all analyses

## Parameter descriptions
# n -- the number of trees to simulate
# seed -- starting number of taxa in random seed tree (must be >= 2)
# stop.by -- aim for number of taxa (n) or amount of time (t) in end tree
# stop.at -- target number of taxa or time of end tree
# leeway -- percentage buffer around target
# birth -- how many births per unit of branch length?
# death -- how many deaths per unit of branch length?
# bias -- what type of ED? 'PE', 'ES' or 'FP'
# min.strength -- min power determing the effect of the bias
# max.strength -- max power determing the effect of the bias

## Shared functions
closeDevices <- function () {
  # make sure all graphical devices are off
  while (!is.null (dev.list ())) {
    dev.off ()
  }
}

## Meta parameter set-up for 4 analyses
n.analyses <- 4 # make sure there are four elements in each meta parameter
meta.n <- rep (2, n.analyses)
meta.seed <- c (2, 2, 2, 100)
meta.birth <- c (2, 2, 2, 1)
meta.death <- c (1, 1, 1, 1)
meta.bias <- c ('FP', 'FP', 'FP', 'FP')
meta.stop.by <- c ('n', 'n', 'n', 't')
meta.stop.at <- c (10, 50, 100, 10)
meta.leeway <- c (10, 10, 10, 10)
meta.min.strength <- c (-1, -1, -1, -1)
meta.max.strength <- c (1, 1, 1, 1)

## If there isn't a results folder, create one
if (!file.exists ('results')) {
  dir.create ('results')
}

## Loop through each analysis and run
for (i in 1:n.analyses) {
  # print analysis number
  cat ('\n--------------------------------')
  cat (paste0 ('\n          Analysis [', i, ']'))
  cat ('\n--------------------------------\n')
  # parameter set-up
  n <- meta.n[i]
  seed <- meta.seed[i]
  birth <- meta.birth[i]
  death <- meta.death[i]
  bias <- meta.bias[i]
  stop.by <- meta.stop.by[i]
  stop.at <- meta.stop.at[i]
  leeway <- meta.leeway[i]
  min.strength <- meta.min.strength[i]
  max.strength <- meta.max.strength[i]
  # results folder set-up
  res.dir <- paste0 ('parameter_set_', i)
  if (!file.exists (file.path ('results', res.dir))) {
    dir.create (file.path ('results', res.dir))
  }
  res.dir <- file.path ('results', res.dir)
  # create a run log
  runlog <- file.path (res.dir, 'runlog.csv')
  if (file.exists (runlog)) {
    file.remove (runlog)
  }
  headers <- data.frame ('treefilename', 'strength', 'bias',
                         'birth', 'death')
  write.table (headers, runlog, sep = ',', row.names = FALSE,
               col.names = FALSE)
  # run
  min.n <- stop.at - (stop.at*leeway/100)
  max.n <- stop.at + (stop.at*leeway/100)
  for (j in 1:n) {
    # print statement
    cat (paste0 ('\n--- Working on model [', j,'] of [', n,
                 '] ---'))
    strength <- runif (1, min.strength, max.strength)
    stop.at <- round (runif (1, min.n, max.n))
    source (file.path ('stages', 'model_trees.R'),
            print.eval = TRUE)
  }
  cat ('\n--- Comparing trees to natural trees ---')
  source (file.path ('stages','compare_trees.R'), print.eval = TRUE)
  cat ('\n------ Model completed ------\n')
}