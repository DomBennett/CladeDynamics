## 17/07/2014
## D.J. Bennett
## Calculating imbalance for real trees

## Libs
library (ape)
library (apTreeshape)
library (caper)
library (geiger)

## Dirs
tree <- read.tree (file.path ('data', 'squamates.tre'))

calcTreeShapeStats <- function (tree) {
  ## Calculates a variety of tree shape statistics and compares
  ##  to stats generated for an equally sized Yule tree
  ## Stats greater than 1 are greater than the Yule reference
  # reference Yule tree
  ref.tree <- sim.bdtree (b = 1, d = 0, n = length (tree$tip.label))
  # convert trees to apTreeshape objects
  ts.tree <- as.treeshape(tree)
  ts.ref.tree <- as.treeshape(ref.tree)
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
  if (is.ultrametric (tree)) {
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
  } else {
    gamma.stat <- tc.stat <- NA
  }
  # return list of stats
  list (colless = colless.stat, sackin = sackin.stat,
        iprime = iprime.stat, gamma = gamma.stat,
        tc = tc.stat)
}