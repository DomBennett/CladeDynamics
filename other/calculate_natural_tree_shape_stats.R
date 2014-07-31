## 31/07/2014
## D.J. Bennett
## Calculate tree shape stats for natural trees

## Parameters
max.n = 1000
min.n = 500

## Libraries
source (file.path ('tools', 'compare_tools.R'))

## Functions
getClades <- function (tree, min.n, max.n) {
  countChildren <- function (node) {
    # add node information to lists children and n
    these.children <- getChildren (tree, node)
    this.n <- length (these.children)
    if (this.n < max.n & this.n > min.n) {
      children <<- c (children, list (these.children))
      n <<- c (n, this.n)
      node.number <<- c (node.number, node)
    }
  }
  checkNode <- function (i, these.names) {
    # return True if names do not overlap
    other.names <- children[[i]]
    !any (these.names %in% other.names)
  }
  # how many nodes to loop through
  ntips <- getSize (tree)
  nodes <- ntips:(ntips + tree$Nnode)
  # create a list for tip names for each clade
  children <- list ()
  # create a vector of n tips and node number for each clade
  n <- node.number <- NULL
  # loop through nodes writing info to children and n
  m_ply (.data = data.frame (node = nodes), .fun = countChildren)
  if (is.null (n)) {
    stop ('No clades found between [', min.n,'] and [', max.n,']')
  }
  # out of those clades, find a non-redundant combination
  # start with a bools matrix that records number of non-redudant
  #  tips between clades
  bools <- matrix (nrow = length (n), ncol = length (n))
  # loop through each node comparing it to others
  for (i in 1:length (n)) {
    # check how many other nodes i node overlaps with
    is <- (1:length (n))[-i]
    these.names <- children[[i]]
    bool <- mdply (.data = data.frame (i = is), .fun = checkNode,
                   these.names)[ ,2]
    # add results to bools matrix
    template <- rep (0, length (n))
    template[i] <- n[i]
    template[is[bool]] <- n[is[bool]]
    # add NAs for combs that have already been recorded
    template[seq (0, i-1)[-1]] <- rep (NA, i-1)
    bools[i, ] <- template
  }
  # sum the rows of the bool matrix
  row.sums <- rowSums (bools, na.rm = TRUE)
  # find the combination with the highest n represented
  best.row <- bools[row.sums == max (row.sums), ]
  # if best.row is nrow () > 2, choose at random
  if (!is.null (dim (best.row))) {
    best.row <- best.row[sample (1:nrow (best.row), 1)]
  }
  best.nodes <- node.number[best.row != 0 & !is.na (best.row)]
  # extract clades and return
  trees <- list ()
  for (best.node in best.nodes) {
    clade.tree <- extract.clade (tree, node = best.node)
    trees <- c (trees, list (clade.tree))
  }
  if (length (trees) > 1) {
    class (trees) <- 'multiPhylo'
    return (trees)
  } else {
    return (trees[[1]])
  }
}

## Dirs
data.dir <- 'data'

## Input
natural.trees <- list ()
tree.files <- list.files (data.dir, '.tre')
for (i in 1:length (tree.files)) {
  tree <- read.tree (file.path (data.dir, tree.files[i]))
  # choose the first tree if list
  if (class (tree) == 'multiPhylo') {
    tree <- tree[[1]]
  }
  # if the tree is bigger than 1000 extract clades that are
  #  greater than 500 and less than 1000
  if (getSize (tree) >= max.n) {
    clade.trees <- getClades (tree, min.n, max.n)
    for (clade.tree in clade.trees) {
      natural.trees <- c (natural.trees, list (clade.tree))
    }
  } else {
    natural.trees <- c (natural.trees, list (tree))
  }
}

## Calculate tree stats
natural.tree.stats <- calcTreeShapeStats (natural.trees)

## Save tree stats
save (natural.tree.stats, file = file.path (data.dir, 'natural_tree_stats.Rd'))