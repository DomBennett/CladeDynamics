## D.J. Bennett
## 21/10/2014
## Clades tools

## Libraries
library (MoreTreeTools)
library (phytools)

.countChildren <- function (tree) {
  # Count the number of extant children for every node
  .count <- function (node.label) {
    node <- which (tree$node.label == node.label)
    node <- length (tree$tip.label) + node
    if (sum (tree$edge[ ,1] == node) > 1) {
      # if the node has two descending edges
      node.children <- getChildren (tree, node)
      extant <- node.children[node.children %in% extants]
      return (length (extant))
    }
    return (0)
  }
  # all internal nodes
  node.labels <- tree$node.label
  # make sure only extant clades are counted
  extants <- getExtant (tree)
  res <- mdply (.data = data.frame (node.label = node.labels),
                .fun = .count)
  colnames (res) <- c ('node', 'n.children')
  res
}

.reformat <- function (clade.performance, sample,
                       progress.bar = 'none') {
  # Take list of list of clade performances and convert
  #  to a dataframe.
  .getTime <- function (i, node) {
    # get success for node at a time point
    data <- clade.performance[[i]]
    if (any (data$node == node)) {
      return (data[data$node == node, 2])
    }
    0
  }
  .getNode <- function (node) {
    mdply (.data = data.frame (
      i = 1:length (clade.performance)),
      .fun = .getTime, node)[ ,2]
  }
  .addNode <- function (node) {
    # add success for node at all time points
    # for a res dataframe
    node <- as.character (node)
    node.success <- .getNode (node)
    res[node] <- node.success
    res <<- res
  }
  # get nodes across times
  nodes <- unique (unlist (llply (.data = clade.performance,
                                  .fun = function (x) as.vector(x$node))))
  # build res dataframe by adding first results
  res <- data.frame (.getNode (nodes[1]))
  colnames (res) <- nodes[1]
  nodes <- data.frame (node = nodes[-1])
  # add to res
  m_ply (.data = nodes, .fun = .addNode, .progress = progress.bar)
  rownames (res) <- seq (sample, sample*nrow (res), sample)
  res
}

getCladeSuccess <- function (trees, sample) {
  ## Calculate the success of clades from a sample of trees
  .get <- function (i) {
    .countChildren (trees[[i]])
  }
  clade.successes <- mlply (.data = data.frame (i = 1:length (trees)),
                            .fun = .get)
  .reformat (clade.successes, sample)
}

calcCladeStats <- function (clades) {
  ## Take dataframe of clade sizes through time, return stats
  ##  for each clade
  .calc <- function (name) {
    clade <- clades[ ,as.character (name)]
    # basic stats
    tot.size <- sum (clade)
    max.size <- max (clade)
    if (tot.size == 0) {
      start <- end <- time.span <- NA
    } else {
      # get start and end
      start <- which (clade > 0)[1]
      end <- which (clade > 0)[sum (clade > 0)]
      # only count non-zero times
      clade <- clade[clade != 0]
      time.span <- length (clade)
    }
    if (time.span < 2 || is.na (time.span)) {
      cm <- cg <- NA
    } else {
      # calculate nominal times
      times <- seq (0, 1, 1/(length (clade) - 1))
      # calculate CM
      cm <- sum (clade * times)/sum (clade)
      # calculate CG
      cg <- sum (clade * (cm - times)^2)/sum (clade)
    }
    data.frame (tot.size, max.size, time.span, start, end, cm, cg)
  }
  mdply (.data = data.frame (name = colnames (clades)),
         .fun = .calc)
}