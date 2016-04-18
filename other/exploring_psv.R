library(MoreTreeTools)
n <- 16
tree <- compute.brlen (stree (n, 'balanced'))
tree$root.edge <- 0.1  # add a root edge
plot (reGravitise (tree, factor = -0.9), show.tip.label=FALSE,
      edge.width=2, root.edge=TRUE)
plot (reGravitise (tree, factor = 0.9), show.tip.label=FALSE,
      edge.width=2, root.edge=TRUE)

low <- reGravitise (tree, factor = -0.9)
high <- reGravitise (tree, factor = 0.9)
calcPSV <- function (tree) {
  samp <- matrix (rep (1, getSize (tree) * 2), nrow=2)
  colnames (samp) <- tree$tip.label
  psv.res <- psv (samp, tree)
  psv.res[1,1]
}
calcPSV(low)
calcPSV(high)


par (mfrow = c (2, 2))
# line tree, high gravity, low PSV
ref.tree <- stree (16, 'left')
ref.tree$edge.length <- rep(0.000001, nrow(ref.tree$edge))
ref.tree$edge.length[c(1,2)] <- 1000
psv.val <- signif(calcPSV(ref.tree), 1)
plot (ref.tree, type='unrooted', show.tip.label = FALSE,
      main=paste0('PSV = ', psv.val))
plot (ref.tree, show.tip.label = FALSE)

# star tree, low gravity, high PSV
tip.edges <- ref.tree$edge[ ,2] %in% 1:length(ref.tree$tip.label)
ref.tree$edge.length[!tip.edges] <- 0.1
ref.tree$edge.length[tip.edges] <- 100
psv.val <- signif(calcPSV(ref.tree), 1)
plot (ref.tree, type='unrooted', show.tip.label = FALSE,
      main=paste0('PSV = ', psv.val))
plot (ref.tree, show.tip.label = FALSE)
