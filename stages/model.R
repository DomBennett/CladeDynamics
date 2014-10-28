## 07/08/2014
## D.J. Bennett
## Model trees using EDBMM

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Parameters
if (!exists ('pars')) {
  print ('here')
  pars <- list (n.model = 2, seed = 2,
                max.birth = 5, min.birth = 1.1,
                max.death = 1, min.death = 1,
                bias = 'FP', stop.by = 'n',
                max.ntaxa = 200, min.ntaxa = 50,
                min.psi = -1, max.psi = 1,
                reference = TRUE, iterations = 100)
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
  headers <- data.frame ('treefilename', 'psi', 'bias',
                         'birth', 'death','ntaxa')
  write.table (headers, runlog, sep = ',',
               row.names = FALSE, col.names = FALSE)
  runlog
}

genParameters <- function (pars) {
  ## Generate parameters for model itreation
  pars$psi = runif (1, pars$min.psi, pars$max.psi)
  pars$birth = runif (1, pars$min.birth, pars$max.birth)
  pars$death = runif (1, pars$min.death, pars$max.death)
  pars$ntaxa = round (runif (1, pars$min.ntaxa, pars$max.ntaxa))
  pars
}

addEntry <- function (runlog, pars) {
  ## Add an entry to runlog.csv
  counter <- nrow (read.csv (runlog))
  treefilename <- paste0 ('tree', counter, '.tre')
  parameters <- data.frame (treefilename,
                            psi = pars$psi,
                            bias = pars$bias,
                            birth = pars$birth,
                            death = pars$death,
                            ntaxa = pars$ntaxa)
  write.table (parameters, runlog, sep = ',', append = TRUE,
               col.names = FALSE, row.names = FALSE)
  treefilename
}

iterateModel <- function (j, pars, runlog, ...) {
  ## Iterate through models
  cat (paste0 ('\n... working on model [', j,'] of [',
               pars$n.model, ']'))
  # set parameters for iteration
  pars <- genParameters (pars)
  treefilename <- addEntry (runlog, pars)
  tree <- runEDBMM (birth = pars$birth,
                    death = pars$death,
                    stop.at = pars$ntaxa,
                    stop.by = pars$stop.by,
                    psi = pars$psi,
                    bias = pars$bias,
                    fossils = FALSE, record = FALSE)
  write.tree (tree, file = file.path (dirname (runlog),
                                      treefilename))
}

## Run
runlog <- dirSetup (name)
m_ply (.data = data.frame (j = 1:pars$n.model),
       .fun = iterateModel, pars, runlog)