## 17/06/2014
## D.J. Bennett
## The rise and fall of clades
## Custom functions for modelling tree growth

## Deps
library (ape)
library (MoreTreeTools)
library (plyr)

## Functions
seedTree <- function (n, age) {
  ## Create a seed tree for growing with n tips
  tree <- rtree (n)
  # add edge lengths
  tree <- compute.brlen (tree)
  tree.age <- getAge (tree, node = length (tree$tip.label) + 1)
  # re-calculate edge lengths based on requested tree age
  tree$edge.length <- tree$edge.length/(tree.age/age)
  # add labels
  tree$node.label <- paste0 ('n', 1:tree$Nnode)
  tree
}

countChildren <- function (tree) {
  # Count the number of extant children for every node
  .count <- function (node.label) {
    node <- which (tree$node.label == node.label)
    node <- length (tree$tip.label) + node
    if (sum (tree$edge[ ,1] == node) > 1) {
      # if the node has two descending edges
      node.children <- getChildren (tree, node)
      extant <- node.children[!node.children %in% extinct]
      return (length (extant))
    }
    return (0)
  }
  # all internal nodes
  node.labels <- tree$node.label
  res <- mdply (.data = data.frame (node.label = node.labels),
                .fun = .count)
  colnames (res) <- c ('node', 'n.children')
  res
}

reformat <- function (clade.performance, interval) {
  # take list of list and convert to a dataframe
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
  m_ply (.data = nodes, .fun = .addNode, .progress = 'time')
  rownames (res) <- seq (interval, interval*nrow (res), interval)
  res
}

calcAddBool <- function (seed.tree, birth, death, max.age, sample.unit = 1) {
  ## Calculate the adding and dropping of tips for
  ##  specified maximum tree age and sample trees produced at set
  ##  units of tree branch growth
  nspp <- length (tree$tip.label)
  tree.age <- getAge (tree, node = length (tree$tip.label) + 1)
  res <- list () # collect a list of addbools for each sampling unit
  while (tree.age < max.age) {
    branch.growth <- 0
    add.bool <- c ()
    while (branch.growth < sample.unit) {
      add.bool <- c (add.bool, sample (c (TRUE, FALSE), size = 1,
                                       prob = c (birth, death)))
      if (add.bool[length (add.bool)]) {
        nspp <- nspp + 1
        branch.growth <- branch.growth + 1/nspp
      }
    }
    tree.age <- tree.age + branch.growth
    res <- c (res, list (add.bool))
  }
  res
}

growTree <- function (add.bool, bias = c ('none', 'PE', 'FP')) {
  # Grow a tree using an equal rates markov model
  #  with specified births and deaths choose species
  #  to speciate based on bias
  randomTipIndex <- function (extant.tips, add = TRUE) {
    # Return a random tip index based on bias and whether
    #  it is being added or not
    if (bias == 'none') {
      return (sample (extant.tips, 1))
    }
    # calculate probs based on bias
    if (bias %in% c ('PE', 'iPE')) {
      # Pendant edge uses all tip edges
      tip.edges <- which (tree$edge[, 2] %in%
                            extant.tips)
      probs <- tree$edge.length[tip.edges]
    } else {
      # Fair proportion uses the proportion of lost branch
      # if the species were lost
      # drop extinct species
      extant.tree <- drop.tip (tree, tip = extinct)
      # calc fp
      probs <- calcFairProportion (extant.tree)
      # reorder to match extant.tips
      probs <-
        probs[match (tree$tip.label[extant.tips], probs[ ,1]), 2]
    }
    # inverse probabilities if adding
    # so that the least evolutionary have the highest chance
    # of speciating (or inverse if iPE or iFP)
    if (bias %in% c ('FP', 'PE')) {
      if (add) {
        probs <- 1/probs
      }
    } else {
      if (!add) {
        probs <- 1/probs
      }
    }
    probs <- probs/sum (probs)
    # return tip index based on probs
    sample (extant.tips, size = 1, prob = probs)
  }
  add <- function (extant.tips) {
    # time passed is 1/length (extant.tips) -- so that
    #  births and deaths are scaled to 1 unit of branch length
    time.passed <- 1/length (extant.tips)
    # find tip edges for all extant tips
    tip.edges <- which (tree$edge[ ,2] %in% extant.tips)
    # extant edges grow by 1/length (extant.tips)
    tree$edge.length[tip.edges] <-
      tree$edge.length[tip.edges] + time.passed
    # find random species to speciate
    to.speciate <- randomTipIndex (extant.tips)
    # find its tip edge
    to.speciate <- which (tree$edge[ ,2] %in% to.speciate)
    # new node must have unique labels
    new.node.label <- paste0 ('n', max.node + 1)
    new.tip.label <- paste0 ('t', max.tip + 1)
    # new node.age is always time.passed -- one time step ago.
    node.age <- time.passed
    # add new tip at random edge
    tree <<- addTip (tree = tree, edge = to.speciate,
                     tip.name = new.tip.label,
                     node.age = node.age,
                     node.label = new.node.label)
    # add 1 to globals to keep names unique!
    max.node <<- max.node + 1
    max.tip <<- max.tip + 1
  }
  drop <- function (extant.tips) {
    # Choose a species to add to the extinct list
    to.drop <- randomTipIndex (extant.tips, add = FALSE)
    # to.drop is an index, find corresponding tip label
    to.drop <- tree$tip.label[to.drop]
    extinct <<- c (extinct, to.drop)
  }
  run <- function (add.bool) {
    # find all extant tips -- their index in tip.label
    extant.tips <- which (!tree$tip.label %in% extinct)
    if (add.bool) {
      add (extant.tips)
    } else if (length (extant.tips) > 2) {
      drop (extant.tips)
    } else {
      # At the beginning of a run this might appear a lot...
      #  try increasing the number of births to deaths
      print ('Drop species -- but too few species in tree. Skipping.')
    }
  }
  m_ply (.data = (add.bool = add.bool), .fun = run)
  tree
}