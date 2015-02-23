# 02/02/2015

source (file.path ('tools','model_tools.R'))

tree <- stree (5, 'left')
tree <- rtree (8)
tree$edge.length <- rep (1, nrow (tree$edge))
plot (tree)
tree.2 <- tree

trees.5 <- list ('unbalanced' = tree.1,
                 'mid' = tree.3,
                 'balanced' = tree.2)

tree <- rtree (6)
tree$edge.length <- rep (1, nrow (tree$edge))
plot (tree)
trees.6 <- c (trees.6, list (tree))
bs.6 <- rep (NA, 5)
n <- 6
exp.sd <- 2*((n*log2(n)) - n + 1)
for (i in 1:5) {
  ds <- colSums(cophenetic.phylo (trees.6[[i]]))
  bs.6[i] <- sqrt(sum ((ds - exp.sd)^2))
}
print (bs.6)
(bs.6 * 3)/exp.sd
round (bs.6)



bstat <- function (tree) {
  n <- getSize (tree)
  tree$edge.length <- rep (1, nrow (tree$edge))
  exp.sd <- 2*((n*log2(n)) - n + 1)
  ds <- colSums(cophenetic.phylo (tree))
  b <- sqrt(sum ((ds - exp.sd)^2))
  b/exp.sd
}

tree <- stree (6, 'left')
tree$edge.length <- rep (1, nrow (tree$edge))
bstat (tree)


n <- 6
bs <- rep (NA, 10)
trees <- list ()
for (i in 1:(length (bs)-1)) {
  tree <- rtree (n)
  trees <- c (trees, list (tree))
  bs[i] <- bstat (tree)
}
tree <- stree (n, 'left')
trees <- c (trees, list (tree))
bs[i+1] <- bstat (tree)
pdf ('bstat.pdf')
trees <- trees[order (bs)]
bs <- bs[order (bs)]
for (i in 1:length (trees)) {
  tree <- trees[[i]]
  tree$edge.length <- rep (1, nrow (tree$edge))
  plot (tree, main = bs[i])
}
dev.off ()


lds <- apply (cophenetic.phylo (tree), 2, max)
rtts <- diag (vcv.phylo (tree))
sum (abs (abs (lds - rtts) - log2 (n)))

exp.rtt <- log2 (n)
exp.ld <- exp.rtt * 2
sum (lds - exp.ld) + sum (rtts - exp.rtt)
