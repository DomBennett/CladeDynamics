n <- 16
par (mfrow = c (2, 2), mar = c (1, 1, 1, 1))
tree <- compute.brlen (stree (n, 'balanced'))
tree$root.edge <- 0.1  # add a root edge
plot (reLoad (tree, factor = -0.9), show.tip.label=FALSE,
      edge.width=2, root.edge=TRUE)
plot (reLoad (tree, factor = 0.9), show.tip.label=FALSE,
      edge.width=2, root.edge=TRUE)

low <- reLoad (tree, factor = -0.9)
high <- reLoad (tree, factor = 0.9)
calcPSV <- function (tree) {
  samp <- matrix (rep (1, getSize (tree) * 2), nrow=2)
  colnames (samp) <- tree$tip.label
  psv.res <- psv (samp, tree)
  psv.res[1,1]
}
calcPSV(ref.tree)

ref.tree <- stree (16)
ref.tree$edge.length <- rep (1, 16)
ref.tree$edge.length[1] <- 0.1
calcPSV(ref.tree)
plot (ref.tree, )
ref.tree$root.edge <- 1
