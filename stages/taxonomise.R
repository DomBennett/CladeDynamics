# 04/06/2015
# D.J. Bennett
# Get taxonomic information for each empirical tree

# LIBS
library (MoreTreeTools)

# FUNCTION
resolve <- function (names) {
  # resolve names using the most commmon taxonomic databases
  datasources <- c (1, 3, 4, 11, 12)
  # Catalogue of life, ITIS, NCBI, GBIF and EOL
  res <- taxaResolve (names, datasource=datasources[1])
  for (d in datasources[-1]) {
    if (all (!is.na (res$search.name))) {
      break
    }
    names <- names[is.na (res$search.name)]
    temp.res <- taxaResolve (names, datasource=d)
    if (all (is.na (res$search.name))) {
      res <- temp.res
    } else {
      res <- rbind (res[!is.na (res$search.name), ], temp.res)
    }
  }
  res
}

getRank <- function (taxa.res, rank, additional=NULL) {
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
  res <- mdply (.data=data.frame(i=1:nrow(taxa.res)), .fun=.do)[ ,2]
  if (all (is.na (res)) & !is.null (additional)) {
    # if no names were resolved at rank, use additional info
    taxa.res <- resolve (additional)
    res <- getRank (taxa.res, rank)
  }
  res
}

getMostCommon <- function (line.rank) {
  line.rank <- names (sort (table (line.rank), TRUE))[1]
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

# COUNT NAMES
all.names <- unlist (tiplabels)
all.names <- unique (all.names)
cat ('\n[', length (all.names), '] unique names.', sep='')

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
  res <- resolve (names)  # resolve names
  if (!all (is.na (res$name.string))) {
    orders <- getRank (res, 'order') # get most common name
    order <- getMostCommon (orders)
    classes <- getRank (res, 'class', order)
    class <- getMostCommon (classes)
    phyla <- getRank (res, 'phylum', c (order, class))
    phylum <- getMostCommon (phyla)
  } else {
    phylum <- class <- order <- NA
  }
  row.element <- data.frame (treefile=treefiles[i], phylum=phylum,
                             class=class, order=order)
  taxoinfo <- rbind (taxoinfo, row.element)
}
taxoinfo <- taxoinfo[-1, ]
cat ('\nDone.')

# OUTPUT
write.csv (taxoinfo, file=file.path (output.dir, 'taxoinfo.csv'),
           row.names=FALSE)