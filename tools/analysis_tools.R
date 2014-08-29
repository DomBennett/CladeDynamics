## 17/06/2014
## D.J. Bennett
## Tools for analysing EDBMM tree runs

## Deps
library (MoreTreeTools)

## Functions
findNonZeros <- function (element) {
  ## Find all results that are not zero
  ##  save for those that occur 1 before
  ##  and 1 after
  bool <- element > 0
  # add shifted bools to get a 1,2 or 3 representing the
  #  position of the number in the sequence
  n.pos <- bool + c (bool[-1], bool[1]) +
    c (bool[length (bool)], bool[-length (bool)])
  # return a bool of all n.pos above 0
  n.pos > 0
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

normalise <- function (element, min.size, min.time) {
  ## Return a list of normalised clade success and time
  ##  centred at its peak
  if (element[length (element)] != 0 |
        element[1] != 0) { # only extinct clades
    return (NA)
  }
  if (max (element) < min.size) {
    return (NA)
  }
  if (sum (element != 0) < min.time) {
    return (NA)
  }
  # remove excessive zeros
  y <- element[findNonZeros (element)]
  # find peak(s)
  peak <- which (y == max (y))
  peak <- round (mean (peak))
  # centre time around peak
  x <- seq ((0 - (peak-1)), (length (y) - peak))
  # normalise
  x <- x/max (x)
  y <- y/max (y)
  list (x = x, y = y)
}

plotNormalisedSuccess <- function (res, min.time.span = 5, min.size = 2) {
  #   if (is.null (dim (res))) {
  #     stop ('Too few suitable clades -- trying re-running')
  #   }
  n.reps <- ceiling (ncol (res)/3)
  cols <- rep (rainbow (3, alpha = 0.7), n.reps)
  plot (x = 0, y = 0, type = 'n', xlim = c (-1,1),
        ylim = c (0, 1), main =
          'Normalised Clade Success Through Time',
        xlab = 'Normalised Time', ylab = 'Normalised Success')
  all.x <- all.y <- c ()
  for (i in 1:ncol (res)) {
    normalised.element <- normalise (res[ ,i], min.size, min.time.span)
    if (any (is.na (normalised.element))) {
      next
    }
    lines (x = normalised.element$x, y = normalised.element$y,
           pch = 19, col = cols[i])
    all.x <- c (all.x, normalised.element$x)
    all.y <- c (all.y, normalised.element$y)
  }
  if (any (is.null (all.x))) {
    cat ('\n... No clades met minimum criteria for normalised plotting')
  } else {
    # add a mean line...
    all.x <- round (all.x, digits = 1)
    uni.all.x <- sort (unique (all.x))
    meaned.y <- rep (NA, length (uni.all.x))
    for (i in 1:length (uni.all.x)) {
      meaned.y[i] <- mean (all.y[all.x %in% uni.all.x[i]])
    }
    lines (x = uni.all.x, y = meaned.y, lwd = 2)
  }
}

.countChildren <- function (tree, extinct = NULL) {
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
    clade <- clades[ ,name]
    # get start and end
    start <- which (clade > 0)[1]
    end <- which (clade > 0)[sum (clade > 0)]
    # reduce clade times to those when it exists
    clade <- clade[clade != 0]
    # basic stats
    tot.size <- sum (clade)
    max.size <- max (clade)
    time.span <- length (clade)
    # calculate nominal times
    times <- seq (0, 1, 1/(length (clade) - 1))
    # calculate CM
    cm <- sum (clade * times)/sum (clade)
    # calculate CG
    cg <- sum (clade * (cm - times)^2)/sum (clade)
    data.frame (tot.size, max.size, time.span, start, end, cm, cg)
  }
  mdply (.data = data.frame (name = colnames (clades)),
                .fun = .calc)
}