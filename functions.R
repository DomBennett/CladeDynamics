## 17/06/2014
## D.J. Bennett
## The clades rise and fall
## Custom functions for modelling tree growth,
##  counting success and plotting success of clades

countChildren <- function (tree) {
  # Count the number of children for every node
  .count <- function (node.label) {
    node <- which (tree$node.label == node.label)
    node <- length (tree$tip.label) + node
    length (getChildren (tree, node))
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

findNonZeros <- function (row) {
  ## Find all results that are not zero
  ##  save for those that occur 1 before
  ##  and 1 after
  non.zero <- row != 0
  first.zero <- which (non.zero)[1]
  if (first.zero > 1) {
    non.zero[first.zero - 1] <- TRUE
  }
  last.zero <- which (non.zero)[sum (non.zero)]
  if (last.zero < length (non.zero)) {
    non.zero[last.zero + 1] <- TRUE
  }
  non.zero
}

plotSuccess <- function (res, section = 'all',
                         thinning = 'all') {
  # Plot as separate lines for the success of each
  # clade through time
  # section = 'all' -- for all clades
  # else do for bottom number
  # thinning allows only a proportion of the
  # clade results to be plotted
  if (section != 'all') {
    top <- ncol (res)
    bottom <- top - section
    res <- res[bottom:top]
  }
  if (thinning != 'all') {
    # choose those to plot by biasing towards the
    # older ones -- makes it easier to see
    # too many young ones
    random.is <- sample (1: ncol (res),
                         thinning, prob = 1/1: ncol (res))
    res <- res[ ,random.is]
  }
  n.clades <- ncol (res)
  # re-use the first 3 rainbow colours
  # red, green and blue
  n.reps <- ceiling (n.clades/3)
  cols <- rep (rainbow (3, alpha = 0.7), n.reps)
  x.max <- max (as.numeric (rownames (res)))
  x.min <- min (as.numeric (rownames (res)))
  y.max <- max (res)
  y.min <- min (res)
  plot (0, type = 'n', xlim = c (x.min, x.max), ylim =
          c (y.min, y.max), ylab = 'N', xlab = 'Time',
        main = 'Clade success through time')
  for (i in 1:n.clades) {
    # don't plot all 0s -- otherwise very cluttered
    # but do plot the zero for start and end of clade
    non.zero <- findNonZeros (res[ ,i])
    y <- res[non.zero, i]
    x <- as.numeric (rownames (res)[non.zero])
    lines (y = y, x = x, col = cols[i], pch = 19)
  }
}

plotNormalisedSuccess <- function (res, min.time.span = 5) {
  ## find all clades that go from 0 to 0
  ## And clades that have minimum time span
  res <- res[ ,colSums (res != 0) > min.time.span]
  res <- res[ ,res[1, ] == 0 & res[nrow (res), ] == 0]
  if (is.null (dim (res))) {
    stop ('Too few suitable clades -- trying re-running')
  }
  n.reps <- ceiling (ncol (res)/3)
  cols <- rep (rainbow (3, alpha = 0.7), n.reps)
  plot (x = 0, y = 0, type = 'n', xlim = c (0.1,1),
        ylim = c (0.1, 1.0), main =
          'Normalised Clade Success Through Time',
        xlab = 'Normalised Time', ylab = 'Normalised Success')
  all.x <- all.y <- c ()
  for (i in 1:ncol (res)) {
    # normlise time steps and success
    non.zero <- findNonZeros (res[ ,i])
    normalised.time.steps <- 1/sum (non.zero)
    x <- seq (normalised.time.steps,
              sum (non.zero)*normalised.time.steps,
              normalised.time.steps)
    y <- res[non.zero,i]/max (res[non.zero,i])
    lines (x = x, y = y, pch = 19, col = cols[i])
    all.x <- c (all.x, x)
    all.y <- c (all.y, y)
  }
  # add a mean line...
  all.x <- round (all.x, digits = 1)
  uni.all.x <- sort (unique (all.x))
  meaned.y <- rep (NA, length (uni.all.x))
  for (i in 1:length (uni.all.x)) {
    meaned.y[i] <- mean (all.y[all.x %in% uni.all.x[i]])
  }
  lines (x = uni.all.x, y = meaned.y, lwd = 2)
}

growTree <- function (iterations, birth = 0.6, death = 0.4) {
  add <- function () {
    # an equal rates markov model:
    # choose a species at random
    # simulate it speciating by adding 1 to all tip edges
    # adding a new tip at a new node in the random species' connecting
    # edge at 1 time step ago
    # choose a species at random
    # grow tip.edges by 1
    tip.edges <- which (tree$edge[, 2] %in%
                          1:length (tree$tip.label))
    tree$edge.length[tip.edges] <-
      tree$edge.length[tip.edges] + 1
    #probs <- tree$edge.length[tip.edges] # weight by branch length?
    random.edge <- sample (x = tip.edges, size = 1)
    new.node.label <- paste0 ('n', max.node + 1)
    new.tip.label <- paste0 ('t', length (tree$tip.label) + 1)
    # new node.age is always 1 -- 1 time step ago.
    node.age <- 1
    # add new tip at random edge
    tree <<- addTip (tree = tree, edge = random.edge,
                     tip.name = new.tip.label,
                     node.age = node.age,
                     node.label = new.node.label)
    max.node <<- max.node + 1
  }
  drop <- function () {
    to.drop <- sample (tree$tip.label, 1)
    tree <<- drop.tip (tree, to.drop)
  }
  run <- function (iteration) {
    add.bool <- sample (c (TRUE, FALSE), size = 1,
                        prob = c (birth, death))
    if (add.bool) {
      add ()
    } else if (length (tree$tip.label) > 2) {
      drop ()
    } else {
      # At the beginning of a run this might appear a lot...
      #  try increasing the number of births to deaths
      print ('Drop species -- but too few species in tree. Skipping.')
    }
  }
  m_ply (.data = (iteration = 1:iterations), .fun = run)
  tree
}