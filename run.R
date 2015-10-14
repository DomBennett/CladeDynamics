## 07/08/2014
## D.J. Bennett
## Run all analyses

## Timestamp
cat (paste0 ('\nrun.R started at [', Sys.time (), ']'))

## Library
library (plyr)

## Load parameters
parfile <- 'pandf.R'
source (file.path ('parameters', parfile))

# if there isn't a results folder, create one
if (!file.exists ('results')) {
  dir.create ('results')
}
# make sure analysis results folder doesn't already exist
for (analysis.name in names (analysis.parameters)) {
  if (file.exists (file.path ('results', analysis.name))) {
    stop (paste0 ('[', analysis.name, '] already exists in results folder'))
  }
}

## Iteration function
iterateAnalyses <- function (i) {
  # analysis parameters
  name <- names (analysis.parameters)[i]
  pars <- analysis.parameters[[i]]
  cat ('######################################\n')
  cat (paste0 ('\n          Analysis [', name, '], [', i, '/', length (analysis.parameters),']'))
  cat ('######################################\n\n')
  cat ('\n--------------------------------')
  cat ('          Model stage ....\n')
  cat ('\n--------------------------------\n')
  source (file.path ('stages', 'model.R'), print.eval = TRUE, local = TRUE)
  cat ('\n--------------------------------')
  cat ('          Calculate stage ....\n')
  cat ('\n--------------------------------\n')
  source (file.path ('stages','calculate.R'), print.eval = TRUE, local = TRUE)
}

## Run
m_ply (.data = data.frame (i = 1:length (analysis.parameters)),
       .fun = iterateAnalyses)

## Timestamp
cat (paste0 ('\nrun.R finished at [', Sys.time (), ']'))