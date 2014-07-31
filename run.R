## 05/07/2014
## D.J. Bennett
## Run all scripts with different parameters

## Parameter descriptions
# seed.n -- how big should the initial tree be?
# time -- how many units of branch length should the tree grow by?
# burning -- how many units of branch length to throw away at start to account for bias created by seed tree?
# sample -- how often should sampling the success of a tree occur?
# birth -- how many births per unit of branch length?
# death -- how many deaths per unit of branch length?
# bias -- what type of ED? 'PE' or 'FP'
# strength -- power determing the effect of the bias 
# min.time.span -- minimum amount of time a clade exists for it to be normalised plotted
# min.size -- minimum maximum size for normalised plotting
# plot.tree.growth -- create .gif of tree growing, default False as this requires installation of ImageMagik

## Shared functions
closeDevices <- function () {
  # make sure all graphical devices are off
  while (!is.null (dev.list ())) {
    dev.off ()
  }
}

## Parameter set-up
seed.n <- 2
time <- 5
burnin <- time*0.1 # 10% of time
sample <- 0.1
birth <- 0.6
death <- 0.4
min.time.span <- 5
min.size <- 5
plot.tree.growth <- FALSE
bias <- 'FP'
# changing parameters
strengths <- rep (c (-1, -0.5, 0, 0.5, 1), each = 5)

## Check if there is a results folder
if (!file.exists ('results')) {
  dir.create ('results')
}

## Create run log
headers <- data.frame ("res.dir", "strength", "bias", "time",
            "sample", "birth", "death", "seed.n",
            "min.time.span", "min.size")
runlog <- file.path (
  'results', paste0 (
    'run_log_', format (Sys.time (), "%H%M_%d%m%y"), '.csv'))
write.table (headers, runlog, sep = ',', row.names = FALSE, col.names = FALSE)

## Run
for (i in 1:length (strengths)) {
  cat (paste0 ('\n\nWorking on model [', i,']\n'))
  strength <- strengths[i]
  source ('1_model.R', print.eval = TRUE)
  source ('2_analysis.R', print.eval = TRUE)
}
source ('3_compare.R', print.eval = TRUE)