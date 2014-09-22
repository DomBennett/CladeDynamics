## 07/08/2014
## D.J. Bennett
## Model trees using EDBMM

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Parameters
# run outside of run.R for testing
if (!exists ('analysis.parameters')) {
  testset <- list (target = 10, leeway = 10,
                   min.strength = -1.5, max.strength = 1,
                   birth = 1, death = 1, bias = 'FP',
                   stop.by = 't', seed = 2, n.model = 5)
  analysis.parameters <- list (testset = testset)
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
  headers <- data.frame ('treefilename',
                         'strength', 'bias',
                         'birth', 'death',
                         'stop.by', 'target')
  write.table (headers, runlog, sep = ',',
               row.names = FALSE, col.names = FALSE)
  runlog
}

addEntry <- function (runlog, mpars) {
  ## Add an entry to runlog.csv
  counter <- nrow (read.csv (runlog))
  treefilename <- paste0 ('tree', counter, '.tre')
  parameters <- data.frame (treefilename,
                            strength = mpars$strength,
                            bias = mpars$bias,
                            birth = mpars$birth,
                            death = mpars$death,
                            stop.by = mpars$stop.by,
                            target = mpars$target)
  write.table (parameters, runlog, sep = ',', append = TRUE,
               col.names = FALSE, row.names = FALSE)
  treefilename
}

iterateModel <- function (j, mpars, runlog, ...) {
  ## Iterate through models
  cat (paste0 ('\n... working on model [', j,'] of [',
               mpars$n.model, ']'))
  # reset varaible parameters
  mpars$strength <- runif (1, mpars$min.strength,
                                      mpars$max.strength)
  mpars$stop.at <- round (runif (1, mpars$min.n,
                                            mpars$max.n))
  treefilename <- addEntry (runlog, mpars)
  # if seed is greater than 2
  if (mpars$seed > 2) {
    # ... create a seed tree of 2-1
    seed.tree <- runEDBMM (birth = 2, death = 1,
                           stop.at = mpars$seed,
                           stop.by = 'n',
                           strength = mpars$strength,
                           bias = mpars$bias,
                           fossils = FALSE, record = FALSE)
    # reset names so no overlap
    seed.tree$tip.label <- paste0 ('t', 1:getSize (seed.tree))
    seed.tree$node.label <- paste0 ('n', 1:seed.tree$Nnode)
    # continue model on ...
    tree <- runEDBMM (birth = mpars$birth,
                      death = mpars$death,
                      stop.at = mpars$stop.at,
                      stop.by = mpars$stop.by,
                      strength = mpars$strength,
                      bias = mpars$bias,
                      fossils = FALSE, record = FALSE)
  } else {
    tree <- runEDBMM (birth = mpars$birth,
                      death = mpars$death,
                      stop.at = mpars$stop.at,
                      stop.by = mpars$stop.by,
                      strength = mpars$strength,
                      bias = mpars$bias,
                      fossils = FALSE, record = FALSE)
  }
  write.tree (tree, file = file.path (dirname (runlog),
                                      treefilename))
}

iterateAnalysis <- function (i) {
  ## Iterate through analyses
  cat ('\n--------------------------------')
  cat (paste0 ('\n          Analysis [', i, ']'))
  cat ('\n--------------------------------\n')
  runlog <- dirSetup (names (analysis.parameters)[i])
  mpars <- analysis.parameters[[i]]
  mpars$min.n <- mpars$target -
    (mpars$target*mpars$leeway/100)
  mpars$max.n <- mpars$target +
    (mpars$target*mpars$leeway/100)
  m_ply (.data = data.frame (j = 1:mpars$n.model),
         .fun = iterateModel, mpars, runlog)
}

## Run
m_ply (.data = data.frame (i = 1:length (analysis.parameters)),
       .fun = iterateAnalysis)