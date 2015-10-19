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
if (!exists ('name')) {
  name <- 'analysis_5'
  runtime <- 2  # how much longer than the original tree
  sample.rate <- 0.01  # sample rate
  ncpus <- 2
  overwrite <- FALSE
}
registerDoMC (ncpus)

## Dirs
input.dir <- file.path ('results', name)
output.dir <- file.path (input.dir, 'clade_results')
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Functions
check <- function (outfile) {
  # check if clade files already exist
  clades.file <- file.path (output.dir, paste0 (outfile, '_clades.csv'))
  cladesi.file <- file.path (output.dir, paste0 (outfile, '_clades_ind.csv'))
  clade.stats.file <- file.path (output.dir, paste0 (outfile, '_clade_stats.csv'))
  cladei.stats.file <- file.path (output.dir, paste0 (outfile, '_clade_ind_stats.csv'))
  trees.file <- file.path (output.dir, paste0 (outfile, '_recorded.rda'))
  res <- c ('cf'=clades.file, 'csf'=clade.stats.file, 'trees'=trees.file,
            'cfi'=cladesi.file, 'csfi'=cladei.stats.file)
  if (overwrite) {
    return (res)
  }
  if (file.exists (clades.file) & file.exists (clade.stats.file) &
        file.exists (trees.file) & file.exists (cladesi.file) &
        file.exists (cladei.stats.file)) {
    return (NA)
  }
  res
}


## Input
runlog <- read.csv (file.path (input.dir, 'runlog.csv'),
                    stringsAsFactors=FALSE)

## Generate clades
cat ('\nGenerating clades for [', name, '] trees ....', sep='')
counter <- foreach (i=1:nrow (runlog)) %dopar% {
  cat ('\n.... working on [', i, '] of [', nrow (runlog), ']', sep='')
  filename <- runlog$treefilename[i]
  # check if its already been run
  outfiles <- check (sub ('\\.tre', '', filename))
  if (is.na (outfiles[1])) {
    return (0)
  }
  # read in tree
  tree <- read.tree (file.path (input.dir, filename))
  # get age and time to run for
  age <- getSize (tree, 'rtt')
  t.stop <- age*runtime
  # sample every sample rate times
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
  clades.ind <- getCladeSuccess (trees, ind=TRUE)
  # save trees as rda to keep nodelabels
  save (trees, file=outfiles['trees'])
  # remove trees
  rm (trees)
  # calc stats
  clade.stats <- calcCladeStats (clades)
  clade.ind.stats <- calcCladeStats (clades.ind)
  # write out
  write.csv (clades, outfiles['cf'])
  write.csv (clade.stats, outfiles['csf'])
  write.csv (clades.ind, outfiles['cfi'])
  write.csv (clade.ind.stats, outfiles['csfi'])
  1
}
cat ('\nDone. Generated clade stats for [', sum (unlist (counter)),
     '] trees.', sep='')

## End
cat (paste0 ('\nStage `clades` finished at [', Sys.time (), ']\n'))