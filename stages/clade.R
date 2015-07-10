## D.J. Bennett
## 21/10/2014
## What do clade rises and falls looks like?

## Start
cat (paste0 ('\nStage `clades` started at [', Sys.time (), ']\n'))

## Libraries (UNIX ONLY)
library (foreach)
library (doMC)
library (MoreTreeTools)

## Parameters
if (!exists ('pars')) {
  name <- 'analysis_5'
  runtime <- 2  # how much longer than the original tree
  sample.rate <- 0.01  # sample rate
  ncpus <- 8
}
registerDoMC (ncpus)
overwrite <- FALSE

## Dirs
wd <- file.path ('results', name)

## Functions
check <- function (outfile) {
  # check if clade files already exist
  clades.file <- file.path (wd, paste0 (outfile, '_clades.csv'))
  clade.stats.file <- file.path (wd, paste0 (outfile, '_clade_stats.csv'))
  res <- c ('cf'=clades.file, 'csf'=clade.stats.file)
  if (overwrite) {
    return (res)
  }
  if (file.exists (clades.file) & file.exists (clade.stats.file)) {
    return (NA)
  }
  res
}


## Input
runlog <- read.csv (file.path (wd, 'runlog.csv'),
                    stringsAsFactors=FALSE)

## Generate clades
cat ('\nGenerating clades for [', name, '] trees ....', sep='')
counter <- foreach (i=1:nrow (runlog)) %dopar% {
  cat ('\n.... working on [', i, '] of [', nrow (runlog), ']', sep='')
  filename <- runlog$treefilename[i]
  # check if its already been run
  outfiles <- check (sub ('\\.tre', '', filename))
  if (is.na (outfiles[1])) {
    next
  }
  # read in tree
  tree <- read.tree (file.path (wd, filename))
  # get age and time to run for
  age <- getSize (tree, 'rtt')
  t.stop <- age*runtime
  # sample every 0.1 times
  sample <- t.stop*sample.rate
  # rename tip labels and node labels to avoid shared names
  tree$tip.label <- paste0 ('t', 1:getSize (tree))
  tree$node.label <- paste0 ('n', 1:tree$Nnode)
  # run forward with b-d 1, record all trees
  trees <- runEDBMM (birth=1,
                     death=1,
                     stop.at=t.stop,
                     stop.by = 't',
                     bias='FP',
                     record=TRUE,
                     fossils = FALSE,
                     eps=runlog$eps[i], # use original eps and sig
                     sig=runlog$sig[i],
                     sample=sample,
                     seed.tree=tree)
  # get clades from trees
  clades <- getCladeSuccess (trees)
  # remove trees
  rm (trees)
  # calc stats
  clade.stats <- calcCladeStats (clades)
  # write out
  write.csv (clades, outfiles['cf'])
  write.csv (clade.stats, outfiles['csf'])
  1
}
cat ('\nDone. Generated clade stats for [', sum (unlist (counter)),
     '] trees.', sep='')

## End
cat (paste0 ('\nStage `clades` finished at [', Sys.time (), ']\n'))