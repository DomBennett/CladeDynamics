## 07/08/2014
## D.J. Bennett
## Run all analyses
## Parameter descriptions
# n.model -- the number of trees to simulate
# seed -- starting number of taxa in random seed tree (must be >= 2)
# stop.by -- aim for number of taxa (n) or amount of time (t) in end tree
# max.ntaxa -- max ntaxa in a tree
# min.ntaxa -- min ntaxa in a tree
# birth -- how many births per unit of branch length?
# death -- how many deaths per unit of branch length?
# bias -- what type of ED? 'PE', 'ES' or 'FP'
# min.sig -- min power determing the effect on speciation
# max.sig -- max power determing the effect on speciation
# min.eps -- min power determing the effect on extinction
# max.eps -- max power determing the effect on extinction
# reference -- normalise shape stats with a Yule reference?
# iterations -- n trees in Yule distribution if reference is True

## Timestamp
cat (paste0 ('\nrun.R started at [', Sys.time (), ']'))

## Library
library (plyr)

## Analysis parameter declarations
pan <- list (n.model = 1000, seed = 2,
            max.birth = 2, min.birth = 2,
            max.death = 1, min.death = 1,
            bias = 'FP', stop.by = 'n',
            max.ntaxa = 500, min.ntaxa = 50,
            min.sig = -1, max.sig = -1,
            min.eps = -1, max.eps = -1,
            reference = TRUE,
            iterations = 100)
pf <- pan  # four separate runs for the different scenarios
pf$min.sig <- 1
pf$max.sig <- 1
pf$min.eps <- -1
pf$max.eps <- -1
de <- pan
de$min.sig <- -1
de$max.sig <- -1
de$min.eps <- 1
de$max.eps <- 1
eph <- pan
eph$min.sig <- 1
eph$max.sig <- 1
eph$min.eps <- 1
eph$max.eps <- 1
analysis.parameters <- list (pan=pan, pf=pf, de=de, eph=eph)
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
# make sure each analysis has reference trees to compare with
data.dir <- file.path ('data', 'treestats')
for (i in 1:length (analysis.parameters)) {
  analysis.name <- names (analysis.parameters)[i]
  min.ntaxa <- analysis.parameters[[i]]$min.ntaxa
  max.ntaxa <- analysis.parameters[[i]]$max.ntaxa
  filename <- paste0 ('min', min.ntaxa, '_max', max.ntaxa, '.Rd')
  if (!file.exists (file.path (data.dir, filename))) {
    stop (paste0 ('No real tree stats have been pre-calculated for [',
                  analysis.name, '], run setup.R with [', min.ntaxa,
                  '] min ntaxa and [', max.ntaxa, '] max ntaxa.'))
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
  cat ('          Compare stage ....\n')
  cat ('\n--------------------------------\n')
  #source (file.path ('stages','compare.R'), print.eval = TRUE, local = TRUE)
}

## Run
m_ply (.data = data.frame (i = 1:length (analysis.parameters)),
       .fun = iterateAnalyses)

## Timestamp
cat (paste0 ('\nrun.R finished at [', Sys.time (), ']'))