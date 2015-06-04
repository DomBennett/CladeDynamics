# 04/06/2015
# D.J. Bennett
# Get taxonomic information for each empirical tree

# LIBS
library (MoreTreeTools)

# FUNCTION
getRank <- function (taxa.res, rank) {
  .do <- function (i) {
    ranks <- ranks[[i]]
    lines <- lines[[i]]
    if (is.na (ranks[1])) {
      res <- NA
    } else if (!rank %in% ranks) {
      res <- NA
    } else {
      res <- lines[ranks==rank]
    }
    res
  }
  ranks <- strsplit (taxa.res$rank, '\\|')
  lines <- strsplit (taxa.res$lineage, '\\|')
  mdply (.data=data.frame(i=1:nrow(taxa.res)), .fun=.do)[ ,2]
}

getMostCommon <- function (line.rank) {
  line.rank <- names (sort (table (line.rank)))[1]
  if (is.null (line.rank)) {
    return (NA)
  }
  line.rank
}

# DIRS
input.dir <- file.path ('data', 'parsed_trees')
output.dir <- file.path ('data', 'treestats')

# INPUT
cat ('\nReading in trees ....')
treefiles <- list.files (input.dir, pattern='.tre')
tiplabels <- list ()
for (i in 1:length (treefiles)) {
  cat ('\n.... [', i, '/', length (treefiles), ']', sep='')
  tree <- read.tree (file=file.path (input.dir, treefiles[i]))
  if (class (tree) == 'multiPhylo') {
    names <- tree[[1]]$tip.label
  } else {
    names <- tree$tip.label
  }
  tiplabels[[treefiles[i]]] <- names
}
cat ('\nDone.')

# RESOLVE NAMES
cat ('\nResolving ....')
taxoinfo <- data.frame (treefile=NA, phylum=NA, class=NA,
                        order=NA)
for (i in 1:length (tiplabels)) {
  cat ('\n.... [', i, '/', length (treefiles), ']', sep='')
  names <- tiplabels[[i]]
  if (length (names) > 100) {
    # only use a subset if many names
    names <- sample (names, 100)
  }
  res <- taxaResolve (names)  # resolve names
  phyla <- getRank (res, 'phylum')
  phylum <- getMostCommon (phyla)  # get most common name
  classes <- getRank (res, 'class')
  class <- getMostCommon (classes)
  orders <- getRank (res, 'order')
  order <- getMostCommon (orders)
  row.element <- data.frame (treefile=treefiles[i], phylum=phylum,
                             class=class, order=order)
  taxoinfo <- rbind (taxoinfo, row.element)
}
cat ('\nDone.')

# OUTPUT
write.csv (taxoinfo, file=file.path (output.dir, 'taxoinfo.csv'),
           row.names=FALSE)