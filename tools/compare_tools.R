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

calcTreeShapeStats <- function (tree, reference = TRUE, n.dist = 100) {
  ## Calculates a variety of tree shape statistics. If reference is TRUE
  ##  compares to stats generated for an equally sized Yule tree
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
      # run multi2di n.dist times, calc stats and take means
      run <- function (i) {
        calcStats (multi2di (tree))
      }
      res <- mdply (.data = data.frame (i = 1:n.dist), .fun = run)[ ,-1]
      return (colMeans (res))
    }
  }
  calcBranchingStats <- function (tree) {
    if (is.ultrametric (tree)) {
      # set tree age to 1
      tree.age <- getAge (tree, node = length (tree$tip.label) + 1)
      tree$edge.length <- tree$edge.length/tree.age
      # gamma stat
      gamma.stat <- gammaStat (tree)
      # total cophenetic distance
      tc.stat <- sum (cophenetic (tree))
    } else {
      tc.stat <- gamma.stat <- NA
    }
    c ('tc.stat' = tc.stat, 'gamma.stat' = tc.stat)
  }
  # calculate topology and branching stats
  stats <- c (calcTopologyStats (tree), calcBranchingStats (tree))
  if (reference) {
    # Pure birth reference
    reference <- sim.bdtree (b = 1, d = 0,
                             n = length (tree$tip.label), stop = 'taxa')
    ref.stats <- c (calcTopologyStats (reference),
                    calcBranchingStats (reference))
    # divide stats by ref stats
    stats <- stats/ref.stats
  }
  as.list (stats)
}

calcMeanTreeShapeStats <- function (trees) {
  tree.stats <- list (colless.stat = NULL, sackin.stat = NULL,
                      iprime.stat = NULL, gamma.stat = NULL,
                      tc.stat = NULL)
  for (i in 1:length (trees)) {
    res <- calcTreeShapeStats (trees[[i]])
    for (each in names (tree.stats)) {
      tree.stats[[each]] <- c (tree.stats[[each]],
                               res[[each]])
    }
  }
  # find means and standard deviations
  for (each in names (tree.stats)) {
    mean.res <- mean (tree.stats[[each]], na.rm = TRUE)
    sd.res <- sd (tree.stats[[each]], na.rm = TRUE)
    res <- list (mean.res, sd.res)
    names (res) <- c (paste0 ('mean.', each), paste0 ('sd.', each))
    tree.stats <- c (tree.stats, res)
  }
  tree.stats
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

histClageAges <- function (metadata) {
  ## Take metadata and plot clade ages in hist for each ED strength
  ages <- strength <- NULL
  for (i in 1:nrow (metadata)) {
    res <- read.csv (file = file.path (
      'results', metadata$res.dir[i], 'clades_through_time.csv'))[ ,-1]
    ages <- c (ages, colSums (res != 0))
    strength <- c (strength, rep (metadata$strength[i], ncol (res)))
  }
  clade.ages <- data.frame (ages = ages, strength = strength)
  hist <- ggplot (clade.ages, aes (x = ages, fill = strength))
  hist + geom_density ()
}