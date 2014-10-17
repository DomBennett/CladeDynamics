## 07/08/2014
## D.J. Bennett
## Run all analyses
## Parameter descriptions
# n.model -- the number of trees to simulate
# seed -- starting number of taxa in random seed tree (must be >= 2)
# stop.by -- aim for number of taxa (n) or amount of time (t) in end tree
# target -- target number of taxa or time of end tree
# leeway -- percentage buffer around target
# birth -- how many births per unit of branch length?
# death -- how many deaths per unit of branch length?
# bias -- what type of ED? 'PE', 'ES' or 'FP'
# min.strength -- min power determing the effect of the bias
# max.strength -- max power determing the effect of the bias

## Timestamp
cat (paste0 ('\nrun.R started at [', Sys.time (), ']'))

## Analysis parameter declarations
analysis.1 <- list (n.model = 100, seed = 2,
                    max.birth = 5, min.birth = 1.1,
                    max.death = 1, min.death = 1,
                    bias = 'FP', stop.by = 'n',
                    max.ntaxa = 200, min.ntaxa = 50,
                    min.psi = -1, max.psi = 1)
analysis.parameters <- list (analysis_1 = analysis.1)
rm (analysis.1)
# if there isn't a results folder, create one
if (!file.exists ('results')) {
  dir.create ('results')
}

## Process
source (file.path ('stages', 'model.R'), print.eval = TRUE)
source (file.path ('stages','compare.R'), print.eval = TRUE)

## Timestamp
cat (paste0 ('\nrun.R finished at [', Sys.time (), ']'))