## 31/07/2014
## D.J. Bennett
## Creating a test for determining tree fractalness

## I aim to create a metric that will capture how symmetrical
##  a tree is both in terms of its toplogy and its branching
## I have decided to use a modified Colless test. Colless test
##  calculates the summed abs difference of descendants from
##  all sister clades. I will do the same but calculate it as
##  the abs difference/summed branch length of all sister clades

## Libraries
library (MoreTreeTools)

## Let's start with a fractal tree: the balanced tree
tree <- stree (64, 'left')
tree <- compute.brlen (tree)
# make sure tree is a standard size
tree$edge.length <- tree$edge.length/sum (tree$edge.length)

## Process
# a little function for calculating the stat for a clade tree
calcStat <- function (tree, node) {
  if (node < getSize (tree)) {
    return (0)
  } else {
    clade.tree <- extract.clade (tree, node)
    return (getSize (clade.tree)/
              getSize (clade.tree, 'pd'))
  }
}
# get the sister clades for all internal nodes
sister.nodes <- getSister (tree)
# remove duplicates sister nodes to prevent double counting
i <- 1
while (TRUE) {
  matches <- match (sister.nodes[ ,1], sister.nodes[ ,2])
  if (all (is.na (matches))) {
    break
  }
  sister.nodes <- sister.nodes[-matches[i], ]
  i <- i + 1
}
# calculate abs diffs
abs.diffs <- rep (NA, nrow (sister.nodes))
for (i in 1:nrow (sister.nodes)) {
  tree1.stat <- calcStat (tree, sister.nodes[i,1])
  tree2.stat <- calcStat (tree, sister.nodes[i,2])
  abs.diffs[i] <- abs (tree1.stat - tree2.stat)
}
res <- sum (abs.diffs)