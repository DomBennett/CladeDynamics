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

.calcTreeShapeStats <- function (tree, reference = TRUE, iterations = 100) {
  ## Hidden workhorse function for calcTreeShapeStats
  calcTopologyStats <- function (tree) {
    calcStats <- function (tree) {
      # first convert tree to treeshape object
      ts.tree <- as.treeshape (tree)
      # summed abs difference between sister clades
      colless.stat <- colless (ts.tree)
      # summed number of ancestors for each tip
      sackin.stat <- sackin (ts.tree)
      # fusco test requires creation of a dataframe
      tree$node.label <- NULL
      tree.data <- data.frame (sp = tree$tip.label, nspp =
                                 rep (1, length (tree$tip.label)))
      tree.data <- comparative.data (phy = tree, dat = tree.data,
                                     names.col = sp)
      tree.fusco.res <- fusco.test (tree.data, rich = nspp,
                                    randomise.Iprime = FALSE)
      iprime.stat <- tree.fusco.res$mean.Iprime
      data.frame (colless.stat = colless.stat, sackin.stat = sackin.stat,
            iprime.stat = iprime.stat)
    }
    if (is.binary.phylo (tree)[1] == TRUE) {
      return (unlist (calcStats (tree)))
    } else {
      # run multi2di iterations times, calc stats and take means
      run <- function (i) {
        calcStats (multi2di (tree))
      }
      res <- mdply (.data = data.frame (i = 1:iterations), .fun = run)[ ,-1]
      return (colMeans (res))
    }
  }
  calcBranchingStats <- function (tree) {
    if (!is.null (tree$edge.length)) {
      if (is.ultrametric (tree)) {
        # set tree age to 1
        tree.age <- getAge (tree, node = length (tree$tip.label) + 1)
        tree$edge.length <- tree$edge.length/tree.age
        # gamma stat
        gamma.stat <- gammaStat (tree)
        # total cophenetic distance
        tc.stat <- sum (cophenetic (tree))
        return (c ('tc.stat' = tc.stat, 'gamma.stat' = gamma.stat))
      }
    }
    c ('tc.stat' = NA, 'gamma.stat' = NA)
  }
  calcReference <- function (n) {
    # Calculate equivalent values for a distribution of Yule trees
    #  and return means
    .calc <- function (i) {
      ## TODO: Yule reference takes a long time, find alternative
      reference <- sim.bdtree (b = 1, d = 0, n = n, stop = 'taxa')
      c (calcTopologyStats (reference), calcBranchingStats (reference))
    }
    res <- mdply (.data = data.frame (i = 1:iterations), .fun = .calc)
    colMeans (res[ ,-1])
  }
  # calculate topology and branching stats
  stats <- c (calcTopologyStats (tree), calcBranchingStats (tree))
  if (reference) {
    ref.stats <- calcReference (getSize (tree))
    # divide stats by ref stats (except for gamma)
    bool <- names (stats) != 'gamma.stat'
    stats[bool] <- stats[bool]/ref.stats[bool]
  }
  as.list (stats)
}

calcTreeShapeStats <- function (tree, reference = TRUE,
                                iterations = 100) {
  ## Calculates a variety of tree shape statistics.
  ##  If reference is TRUE compares to stats generated
  ##  for an equally sized Yule tree.
  if (class (tree) == 'multiPhylo' | class (tree) == 'list') {
    stats <- list (colless.stat = NULL, sackin.stat = NULL,
                   iprime.stat = NULL, gamma.stat = NULL,
                   tc.stat = NULL)
    eachTree <- function (i) {
      res <- .calcTreeShapeStats (tree[[i]], reference = reference,
                                  iterations = iterations)
      for (each in names (stats)) {
        stats[[each]] <<- c (stats[[each]],
                             res[[each]])
      }
    }
    m_ply (.data = data.frame (i = 1:length (tree)),
           .fun = eachTree)
    # find means and standard deviations
    for (each in names (stats)) {
      mean.res <- mean (stats[[each]], na.rm = TRUE)
      sd.res <- sd (stats[[each]], na.rm = TRUE)
      res <- list (mean.res, sd.res)
      names (res) <- c (paste0 ('mean.', each),
                        paste0 ('sd.', each))
      stats <- c (stats, res)
    }
    return (stats)
  } else {
    return (.calcTreeShapeStats (tree, reference = reference,
                                 iterations = iterations))
  }
}

extractStat <- function (simulated.tree.stats, stat.name) {
  ## Extract a stat from the simulated.tree.stats list
  ## Return a vector/list of stats in order that they appear in
  ##  list.
  if (grepl ('(mean|sd)', stat.name)) {
    res <- rep (NA, length (simulated.tree.stats))
    for (i in 1:length (simulated.tree.stats)) {
      res[i] <- simulated.tree.stats[[i]][[stat.name]]
    }
  } else {
    res <- list ()
    for (i in 1:length (simulated.tree.stats)) {
      res <- c (res, list (simulated.tree.stats[[i]][[stat.name]]))
    }
  }
  res
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