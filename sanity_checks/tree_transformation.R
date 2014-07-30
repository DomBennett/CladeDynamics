## 25/07/2014
## D.J. Bennett
## Check how tree transformations of branch lengths affects tree stats

## Libs
source (file.path ('tools', 'compare_tools.R'))

## Create different sizes of trees using the same process -- birth-death
## ... along with reference trees. They should all have the same tree shape
## stats, therefore the optimal transformation should produce stats with the least variance
how.many.trees <- 100
ns <- runif (how.many.trees, min = 100, max = 1000)
trees <- list ()
references <- list ()
modelTrees <- function (i) {
  tree <- drop.extinct (sim.bdtree (b = 2, d = 1, n = ns[i],
                                    stop = 'taxa', extinct = FALSE))
  reference <- sim.bdtree (b = 1, d = 0, n = ns[i], stop = 'taxa')
  trees <<- c (trees, list (tree))
  references <<- c (references, list (reference))
}
m_ply (.data = data.frame (i = 1:how.many.trees),
       .fun = modelTrees, .progress = 'time')

## Calculate mean SD for each statistic using different transformations
stat.names <- c ('colless.stat', 'sackin.stat', 'iprime.stat',
                 'gamma.stat', 'tc.stat')
sds.1 <- sds.2 <- sds.3 <- rep (NA, length (stat.names))
names (sds.1) <- names (sds.2) <- names (sds.3) <- stat.names
# Untransformed
trees.stats <- calcTreeShapeStats (trees, reference = FALSE)
refs.stats <- calcTreeShapeStats (references, reference = FALSE)
for (i in 1:length (stat.names)) {
  sds.1[i] <- sd (trees.stats[[stat.names[i]]]/
                    refs.stats[[stat.names[i]]])
}
# Scaled so all tree branches sum to 1
trans.trees.1 <- trans.refs.1 <- list ()
for (i in 1:length (trees)) {
  trans.tree.1 <- trees[[i]]
  trans.tree.1$edge.length <- trans.tree.1$edge.length /
    sum (trans.tree.1$edge.length)
  trans.ref.1 <- references[[i]]
  trans.ref.1$edge.length <- trans.ref.1$edge.length /
    sum (trans.ref.1$edge.length)
  trans.trees.1 <- c (trans.trees.1, list (trans.tree.1))
  trans.refs.1 <- c (trans.refs.1, list (trans.ref.1))
}
trees.stats <- calcTreeShapeStats (trans.trees.1, reference = FALSE)
refs.stats <- calcTreeShapeStats (trans.refs.1, reference = FALSE)
for (i in 1:length (stat.names)) {
  sds.2[i] <- sd (trees.stats[[stat.names[i]]]/
                    refs.stats[[stat.names[i]]])
}
# Scaled so trees have same age
trans.trees.2 <- trans.refs.2 <- list ()
for (i in 1:length (trees)) {
  trans.tree.2 <- trees[[i]]
  trans.tree.2$edge.length <- trans.tree.2$edge.length /
    getAge (trans.tree.2, node = length (trans.tree.2$tip.label) + 1)
  trans.ref.2 <- references[[i]]
  trans.ref.2$edge.length <- trans.ref.2$edge.length /
    getAge (trans.ref.2, node = length (trans.ref.2$tip.label) + 1)
  trans.trees.2 <- c (trans.trees.2, list (trans.tree.2))
  trans.refs.2 <- c (trans.refs.2, list (trans.ref.2))
}
trees.stats <- calcTreeShapeStats (trans.trees.2, reference = FALSE)
refs.stats <- calcTreeShapeStats (trans.refs.2, reference = FALSE)
for (i in 1:length (stat.names)) {
  sds.3[i] <- sd (trees.stats[[stat.names[i]]]/
                    refs.stats[[stat.names[i]]])
}

# Topology stats should be unaffected
sds.1
sds.2
sds.3
# Hmmm I would have thought that the branching stats would be affected.
#  Looks like any transformation is overkill.