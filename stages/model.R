## 07/08/2014
## D.J. Bennett
## Model trees using EDBMM

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Parameters
# run outside of run.R for testing
if (!exists ('analysis.parameters')) {
  testset <- list (n.model = 2, seed = 2,
                   max.birth = 5, min.birth = 1.1,
                   max.death = 1, min.death = 1,
                   bias = 'FP', stop.by = 'n',
                   max.ntaxa = 200, min.ntaxa = 50,
                   min.psi = -1, max.psi = 1)
  analysis.parameters <- list (testset = testset)
  rm (testset)
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
  headers <- data.frame ('treefilename', 'psi', 'bias',
                         'birth', 'death','ntaxa')
  write.table (headers, runlog, sep = ',',
               row.names = FALSE, col.names = FALSE)
  runlog
}

genParameters <- function (mpars) {
  ## Generate parameters for model itreation
  mpars$psi = runif (1, mpars$min.psi, mpars$max.psi)
  mpars$birth = runif (1, mpars$min.birth, mpars$max.birth)
  mpars$death = runif (1, mpars$min.death, mpars$max.death)
  mpars$ntaxa = round (runif (1, mpars$min.ntaxa, mpars$max.ntaxa))
  mpars
}

addEntry <- function (runlog, mpars) {
  ## Add an entry to runlog.csv
  counter <- nrow (read.csv (runlog))
  treefilename <- paste0 ('tree', counter, '.tre')
  parameters <- data.frame (treefilename,
                            psi = mpars$psi,
                            bias = mpars$bias,
                            birth = mpars$birth,
                            death = mpars$death,
                            ntaxa = mpars$ntaxa)
  write.table (parameters, runlog, sep = ',', append = TRUE,
               col.names = FALSE, row.names = FALSE)
  treefilename
}

iterateModel <- function (j, mpars, runlog, ...) {
  ## Iterate through models
  cat (paste0 ('\n... working on model [', j,'] of [',
               mpars$n.model, ']'))
  # set parameters for iteration
  mpars <- genParameters (mpars)
  treefilename <- addEntry (runlog, mpars)
  tree <- runEDBMM (birth = mpars$birth,
                    death = mpars$death,
                    stop.at = mpars$ntaxa,
                    stop.by = mpars$stop.by,
                    psi = mpars$psi,
                    bias = mpars$bias,
                    fossils = FALSE, record = FALSE)
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
  m_ply (.data = data.frame (j = 1:mpars$n.model),
         .fun = iterateModel, mpars, runlog)
}

## Run
m_ply (.data = data.frame (i = 1:length (analysis.parameters)),
       .fun = iterateAnalysis)