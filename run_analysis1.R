## 07/08/2014
## D.J. Bennett
## Analysis 1: How does an ED bias affect tree balance?


## Parameter descriptions
# n.trees -- the number of trees to simulate
# target -- target number of taxa in end tree
# leeway -- percentage buffer around target
# size.classes -- range of tree sizes
# birth -- how many births per unit of branch length?
# death -- how many deaths per unit of branch length?
# bias -- what type of ED? 'PE' or 'FP'
# strength -- power determing the effect of the bias
# min.strength
# max.strength

## Shared functions
closeDevices <- function () {
  # make sure all graphical devices are off
  while (!is.null (dev.list ())) {
    dev.off ()
  }
}

## Parameter set-up
n <- 10
birth <- 2
death <- 1
bias <- 'FP'
target <- 100
leeway <- 10
min.strength <- -1
max.strength <- 1

## Check if there is a results folder
if (!file.exists ('results')) {
  dir.create ('results')
}
if (!file.exists (file.path ('results', 'analysis_1'))) {
  dir.create (file.path ('results', 'analysis_1'))
}
res.dir <- file.path ('results', 'analysis_1')

## Create run log
runlog <- file.path (res.dir, 'runlog.csv')
if (!file.exists (runlog)) {
  file.remove (runlog)
}
headers <- data.frame ('treefilename', 'strength', 'bias',
                       'birth', 'death')
write.table (headers, runlog, sep = ',', row.names = FALSE,
             col.names = FALSE)

## Run
min.n <- target - (target*leeway/100)
max.n <- target + (target*leeway/100)
for (i in 1:n) {
  cat (paste0 ('\n------ Working on model [', i,'] of [', n,
               '] ------'))
  strength <- runif (1, min.strength, max.strength)
  n.taxa <- round (runif (1, min.n, max.n))
  source (file.path ('stages', '1_model_taxa.R'),
          print.eval = TRUE)
}
cat ('\n------ Comparing trees to natural trees ------')
source (file.path ('stages','2_compare.R'), print.eval = TRUE)
cat ('\n------ Completed ------')