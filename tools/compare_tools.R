## 17/07/2014
## D.J. Bennett
## Tools for comparing MRMM tree runs

## Libs
library (MoreTreeTools)
library (ape)
library (apTreeshape)
library (caper)
library (geiger)
library (ggplot2)

getAge <- function (tree, node = 'all') {
  run <- function (node) {
    # function for calculating age of node
    term.node <- length (tree$tip.label) + 1
    # if it's the root node, return tree.age
    if (term.node == node) {
      return (tree.age)
    }
    # if it's a tip return 0
    if (node < length (tree$tip.label) & is.ultrametric (tree)) {
      return (0)
    }
    # else find all its edge.lengths and subtract from tree.age
    edges <- c ()
    while (node != term.node) {
      edges <- c (edges, which (tree$edge[ ,2] == node))
      node <- tree$edge[tree$edge[ ,2] == node, 1]
    }
    return (tree.age - sum (tree$edge.length[edges]))
  }
  tree.age <- max (diag (vcv.phylo (tree)))
  if (node != 'all') {
    return (run (node))
  }
  # else node == all, run on all nodes
  nodes <- 1:(length (tree$tip.label) + tree$Nnode)
  res <- mdply (.data = data.frame (node = nodes), .fun = run)
  colnames (res)[2] <- 'age'
  res
}

calcTreeStats <- function (trees) {
  engine <- function (i) {
    tree <- trees[[i]]
    # imbalance stats
    ts.tree <- as.treeshape (tree)
    colless.stat <- colless (ts.tree, 'yule')
    sackin.stat <- sackin (ts.tree, 'yule')
    # branching stats
    if (is.null (tree$edge.length)) {
      gamma.stat <- tci.stat <- NA
    } else {
      gamma.stat <- gammaStat (tree)
      # weight TCI by total branch length and number of tips
      tci.stat <- sum (cophenetic (tree))/
        sum (tree$edge.length)/length(tree$tip.label)
    }
    data.frame (colless = colless.stat, sackin = sackin.stat,
                gamma = gamma.stat, tci = tci.stat)
  }
  mdply (.data = data.frame (i = 1:length (trees)), .fun = engine)
}

drawCorresPoints <- function (model, distribution) {
  ## Plot mean, mean and 5% and 95% quantiles of an x distribution's
  ##  equivalent y values
  .draw <- function (actual, predicted, ...) {
    lines (x = c (actual, actual, min.x),
           y = c (min.y, predicted, predicted), ...)
  }
  # find the minimum values of x and y in plotted space
  min.y <- floor (min (model$model[ ,'y'])) + floor (min (model$model[ ,'y']))
  min.x <- floor (min (model$model[ ,'x'])) + floor (min (model$model[ ,'x']))
  # get actual
  q1.x <- quantile (distribution, 0.05, na.rm = TRUE)
  q2.x <- quantile (distribution, 0.95, na.rm = TRUE)
  m.x <- mean (distribution, na.rm = TRUE)
  # get predicted values
  q1.y <- predict (model, data.frame (x = q1.x))
  q2.y <- predict (model, data.frame (x = q2.x))
  m.y <- predict (model, data.frame (x = m.x))
  # draw lines
  .draw (q1.x, q1.y, col = 'red', lty = 2)
  .draw (q2.x, q2.y, col = 'red', lty = 2)
  .draw (m.x, m.y, lwd = 2)
}

calcDistDiff <- function (dist.1, dist.2, n = 1000) {
  # Randomly pull from dists, calc difference, calc mean
  .calc <- function (i) {
    sample.1 <- sample (dist.1, 1)
    sample.2 <- sample (dist.2, 1)
    sample.1 - sample.2
  }
  res <- mdply (.data = data.frame (i = 1:n),
                .fun = .calc)[ ,2]
  mean (res)
}

windowAnalysis <- function (sim.stats, real.stats, size,
                            psi, incr = 0.01,
                            scale = TRUE) {
  .calc <- function (i) {
    samp <- sim.stats[psi < maxs[i] &
                        psi > mins[i]]
    data.frame (diff = abs (calcDistDiff (samp, real.stats)),
                mid = (maxs[i]+mins[i])/2)
  }
  maxs <- seq (from = (min (psi) + size), to = max (psi), by = incr)
  mins <- seq (from = min (psi), to = (max (psi) - size), by = incr)
  res <- mdply (.data = data.frame (i = 1:length (maxs)),
                .fun = .calc)[ ,-1]
  if (scale) {
    res[['diff']] <- res[['diff']] / max (res[['diff']])
  }
  res
}