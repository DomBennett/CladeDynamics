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

calcTreeShapeStats <- function (tree) {
  ## Calculates a variety of tree shape statistics and compares
  ##  to stats generated for an equally sized Yule tree
  ## Stats greater than 1 are greater than the Yule reference
  calcImabalanceStats <- function (tree, ref.tree) {
    # first convert trees to treeshape objects
    ts.tree <- as.treeshape (tree)
    ts.ref.tree <- as.treeshape (ref.tree)
    # calculate toplogical measures of imbalance
    # summed abs difference between sister clades
    colless.stat <- colless (ts.tree)/colless (ts.ref.tree)
    # summed number of ancestors for each tip
    sackin.stat <- sackin (ts.tree)/sackin (ts.ref.tree)
    # fusco test requires creation of a dataframe
    tree$node.label <- NULL
    tree.data <- data.frame (sp = tree$tip.label, nspp =
                               rep (1, length (tree$tip.label)))
    tree.data <- comparative.data (phy = tree, dat = tree.data,
                                   names.col = sp)
    tree.fusco.res <- fusco.test (tree.data, rich = nspp,
                                  randomise.Iprime = FALSE)
    ref.tree$node.label <- NULL
    ref.tree.data <- data.frame (sp = ref.tree$tip.label, nspp =
                                   rep (1, length (ref.tree$tip.label)))
    ref.tree.data <- comparative.data (phy = ref.tree, dat = ref.tree.data,
                                       names.col = sp)
    ref.tree.fusco.res <- fusco.test (ref.tree.data, rich = nspp,
                                      randomise.Iprime = FALSE)
    iprime.stat <- tree.fusco.res$mean.Iprime/
      ref.tree.fusco.res$mean.Iprime
    list (colless.stat = colless.stat, sackin.stat = sackin.stat,
          iprime.stat = iprime.stat)
  }
  calcShapeStats <- function (tree, ref.tree) {
    # calculate tree branch based metrics
    # set tree ages to 1
    ref.tree.age <- getAge (ref.tree, node =
                              length (ref.tree$tip.label) + 1)
    ref.tree$edge.length <- ref.tree$edge.length/ref.tree.age
    tree.age <- getAge (tree, node =
                          length (tree$tip.label) + 1)
    tree$edge.length <- tree$edge.length/tree.age
    # measure of point of diversification
    # measure difference for gamma stat
    gamma.stat <- gammaStat (tree) - gammaStat (ref.tree)
    # total cophenetic distance
    tc.stat <- sum (cophenetic (tree))/sum (
      cophenetic (ref.tree))
    list (tc.stat = tc.stat, gamma.stat = gamma.stat)
  }
  # create reference Yule tree
  ref.tree <- sim.bdtree (b = 1, d = 0, n = length (tree$tip.label))
  # calculate imbalance metrics
  if (is.binary.phylo (tree)[1] == TRUE) {
    imbalance.stats <- calcImabalanceStats (tree, ref.tree)
  } else {
    # If tree is not binary create ...
    imbalance.stats <- list (colless.stat = c (), sackin.stat = c (),
                             iprime.stat = c ())
    for (i in 1:100) {
      # ... a random distribution of 100 trees ...
      new.tree <- multi2di (tree)
      temp.imbalance.stats <- calcImabalanceStats (new.tree, ref.tree)
      imbalance.stats$colless.stat <- c (imbalance.stats$colless.stat,
                                         temp.imbalance.stats$colless.stat)
      imbalance.stats$sackin.stat <- c (imbalance.stats$sackin.stat,
                                         temp.imbalance.stats$sackin.stat)
      imbalance.stats$iprime.stat <- c (imbalance.stats$iprime.stat,
                                         temp.imbalance.stats$iprime.stat)
    }
    # ... and take the mean stats of these trees
    imbalance.stats$colless.stat <- mean (imbalance.stats$colless.stat)
    imbalance.stats$sackin.stat <- mean (imbalance.stats$sackin.stat)
    imbalance.stats$iprime.stat <- mean (imbalance.stats$iprime.stat)
  }
  # calculate shape statistics
  if (is.ultrametric (tree)) {
    # if tree is ultrametric, calculate stats
    shape.stats <- calcShapeStats (tree, ref.tree)
  } else {
    # else make NA
    shape.stats <- list (gamma.stat = NA, tc.stat = NA)
  }
  # return list of stats
  c (imbalance.stats, shape.stats)
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