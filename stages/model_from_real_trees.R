## 07/08/2014
## D.J. Bennett
## Model trees using EDBMM from real trees

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Parameters
if (!exists ('pars')) {
  pars <- list (n.model = 10,
                bias = 'FP', stop.by = 't',
                min.sig = -1, max.sig = 1,
                min.eps = -1, max.eps = 1)
  name <- 'testset'
}

## Functions
dirSetup <- function (analysis.name) {
  ## Create a results folder and runlog.csv
  if (!file.exists (file.path ('results', analysis.name))) {
    dir.create (file.path ('results', analysis.name))
  }
  res.dir <- file.path ('results', analysis.name)
  runlog <- file.path (res.dir, 'runlog.csv')
  if (file.exists (runlog)) {
    file.remove (runlog)
  }
  headers <- data.frame ('treefilename', 'reference', 'sig',
                         'eps', 'bias', 'birth', 'death','time')
  write.table (headers, runlog, sep = ',',
               row.names = FALSE, col.names = FALSE)
  runlog
}

addEntry <- function (runlog, pars) {
  ## Add an entry to runlog.csv
  counter <- nrow (read.csv (runlog))
  treefilename <- paste0 ('tree', counter, '.tre')
  parameters <- data.frame (treefilename,
                            reference = pars$real.tree,
                            sig = pars$sig,
                            eps = pars$eps,
                            bias = pars$bias,
                            birth = pars$birth,
                            death = pars$death,
                            time = pars$time)
  write.table (parameters, runlog, sep = ',', append = TRUE,
               col.names = FALSE, row.names = FALSE)
  treefilename
}

readInTrees <- function (treedir) {
  # read in real trees and return all those that are ultrametric
  treeinfo <- read.csv (file.path ('data', 'parsed_trees',
                                   'treeinfo.csv'),
                        stringsAsFactors=FALSE)
  treeinfo <- treeinfo[treeinfo$ultra | treeinfo$chronos, ]
  treeinfo <- treeinfo[1:10,]
  trees <- list ()
  for (i in 1:nrow(treeinfo)) {
    filename <- treeinfo$filename[i]
    tree <- read.tree (file.path (treedir, filename))[[1]]
    # rescale branch lengths to sum to 1
    tree$edge.length <- tree$edge.length/sum (tree$edge.length)
    trees[filename] <- list (tree)
  }
  return(trees)
}

iterateModel <- function (j, pars, runlog) {
  # Iterate through models
  cat (paste0 ('\n...... working on model [', j,'] of [',
               pars$n.model, ']'))
  # set parameters for iteration
  pars$sig = runif (1, pars$min.sig, pars$max.sig)
  pars$eps = runif (1, pars$min.eps, pars$max.eps)
  tree <- runEDBMM (birth = pars$birth,
                    death = pars$death,
                    stop.at = pars$time,
                    stop.by = 't',
                    seed.tree = pars$seed.tree,
                    sig = pars$sig,
                    eps = pars$eps,
                    bias = pars$bias,
                    fossils = FALSE, record = FALSE)
  # record model run in runlog
  treefilename <- addEntry (runlog, pars)
  write.tree (tree, file = file.path (dirname (runlog),
                                      treefilename))
}

iterateRealTrees <- function(i, pars, runlog) {
  cat (paste0 ('\n... working on real tree [', i,'] of [',
               length (real.trees), ']'))
  tree <- real.trees[[i]]
  # set birth and death to 1
  pars$birth <- 1
  pars$death <- 1
  # run for twice the amount of time
  pars$time <- getSize (tree, 'rtt') * 2
  # get real tree name and add to pars
  pars$real.tree <- names (real.trees)[[i]]
  # add tree to pars
  pars$seed.tree <- tree
  # iterate model
  m_ply (.data = data.frame (j = 1:pars$n.model),
         .fun = iterateModel, pars, runlog)
}

## Run
treedir <- file.path ('data', 'parsed_trees')
real.trees <- readInTrees (treedir)
runlog <- dirSetup (name)
m_ply (.data = data.frame (i = 1:length (real.trees)),
       .fun = iterateRealTrees, pars, runlog)

for (i in 1:100){
  tree <- runEDBMM (birth=1, death=1, stop.at=100, eps=0, sig=0)
  print(getSize (tree, 'rtt'))
}
