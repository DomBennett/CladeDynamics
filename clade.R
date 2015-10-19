## 07/08/2014
## D.J. Bennett
## Run multiple clade analyses

## Timestamp
cat (paste0 ('\nclade_run.R started at [', Sys.time (), ']'))

## Library
library (plyr)

## Load parameters
analysis.names <- c ('pan', 'de', 'hydra')  # vector of analysis names to run from
runtime <- 2  # how much longer than the original tree
sample.rate <- 0.01  # sample rate
ncpus <- 8
overwrite <- FALSE

# make sure analysis results folders exist
for (analysis.name in analysis.names) {
  if (!file.exists (file.path ('results', analysis.name))) {
    stop (paste0 ('[', analysis.name, '] doesn\'t exists in results folder'))
  }
}

## Loop
for (name in analysis.names) {
  cat ('######################################\n')
  cat (paste0 ('\n          Clade [', name, ']'))
  cat ('######################################\n\n')
  source (file.path ('stages', 'clade.R'))
}

## Timestamp
cat (paste0 ('\nclade_run.R finished at [', Sys.time (), ']'))