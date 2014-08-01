## 31/07/2014
## D.J. Bennett
## Creating a test for determining tree symmetry

## I aim to create a metric that will capture how symmetrical
##  a tree is both in terms of its toplogy and its branching
## I have decided to use a modified Colless test. Colless test
##  calculates the summed abs difference of descendants from
##  all sister clades. I will do the same but calculate it as
##  the abs difference/summed branch length of all sister clades

## The Function
testSymmetry <- function (tree, reference = FALSE, iterations = 999,
                          return.clade.diffs = FALSE) {
  ## An extension to Colless test. All sister clades are compared
  ##  as: sum(|(Nr/PDr) - (Nl/PDl)|), where N is the number of taxa, PD the total
  ##  branch length and l and r are all the left and right sister clades
  ##  in the tree.
  ## If reference is True, return a dataframe comparing the stat of tree given
  ##  to the same stat calculated for a distribution of Yule trees of the same
  ##  size. P value returned is the proportion of Yule trees that had a stat
  ##  smaller than the tree given -- a measure of how significantly symmetrical
  ##  the tree is. N.B. small trees that are symmetrical may not be significantly
  ##  symmetrical because smaller trees have a small number of possible trees.
  getNRSisters <- function (tree) {
    # Return all sister clades for all internal nodes, without redundancy
    sister.nodes <- getSister (tree)
    # remove duplicates sister nodes to prevent double counting
    i <- 1
    while (TRUE) {
      matches <- match (sister.nodes[ ,1], sister.nodes[ ,2])
      if (all (is.na (matches))) {
        break
      }
      next.to.drop <- which (!is.na (matches))[1]
      sister.nodes <- sister.nodes[-matches[next.to.drop], ]
      i <- i + 1
    }
    sister.nodes
  }
  calcDiff <- function (tree, node) {
    # Return N/PD or zero if clade is of size 1
    if (node <= getSize (tree)) {
      return (0)
    } else {
      clade.tree <- extract.clade (tree, node)
      return (getSize (clade.tree)/
                getSize (clade.tree, 'pd'))
    }
  }
  calcDiffs <- function (tree) {
    .calc <- function (i) {
      tree1.stat <- calcDiff (tree, sister.nodes[i,1])
      tree2.stat <- calcDiff (tree, sister.nodes[i,2])
      abs (tree1.stat - tree2.stat)
    }
    # get NR sisters
    sister.nodes <- getNRSisters (tree)
    # calculate abs diffs
    mdply (.data = data.frame (i = 1:nrow (sister.nodes)),
           .fun = .calc)
  }
  calcMeanRef <- function (n, iterations) {
    # Return a mean stat from a distribution yule trees
    #  of equivalent size
    .calc <- function (i) {
      reference.tree <- sim.bdtree (b = 1, d = 0, stop = 'taxa',
                                    n = n)
      reference.tree$edge.length <- reference.tree$edge.length/
        sum (reference.tree$edge.length)
      sum (calcDiffs (reference.tree))
    }
    ref.diffs <- mdply (.data = data.frame (i = 1:iterations),
                        .fun = .calc)[ ,2]
    ref.diffs
  }
  # Make tree branch lengths == 1
  tree$edge.length <- tree$edge.length/sum (tree$edge.length)
  tree.diffs <- calcDiffs (tree)
  tree.summed.diff <- sum (tree.diffs)
  if (reference) {
    ref.mean.summed.diffs <- calcMeanRef (getSize (tree), iterations)
    ref.grand.mean <- mean (ref.mean.summed.diffs)
    ref.grand.sd <- sd (ref.mean.summed.diffs)
    # stat defined as the tree/ref
    stat <- tree.summed.diff/ref.grand.mean
    # p.value -- is tree significantly more symmetrical than ref?
    #  i.e. what proportion of the refs had a lower summed difference
    p.value <- sum (ref.mean.summed.diffs < tree.summed.diff)/iterations
    # return data.frame
    return (data.frame (stat, tree = tree.summed.diff,
                        ref.mean = ref.grand.mean, ref.sd = ref.grand.sd,
                        p.value))
  } else if (return.clade.diffs) {
    return (tree.diffs)
  } else {
    return (tree.summed.diff)
  }
}

## Libraries
library (MoreTreeTools)
library (geiger) # TODO import sim.bdtree for MTT

## Try out different trees and see what their stats are
tree <- stree (64, 'left')
tree <- compute.brlen (tree)
testSymmetry (tree, reference = TRUE)