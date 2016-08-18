# MAKE DUMMY TREES FOR SHOWING EXPECTED DIFFERENCES IN SHAPE BETWEEN
# DIFFERENT SIMULATION SCENARIOS.

reGravitise <- function(tree, factor) {
  # Rehashed MTT function with correction to prevent negative branch lengths
  calcUp <- function(edge, change, preceding=TRUE) {
    if(preceding) {
      res <- edge + change
    } else {
      res <- edge - change
    }
    res
  }
  calcDwn <- function(edge, change, preceding=TRUE) {
    if(preceding) {
      res <- edge - change
    } else {
      res <- edge + change
    }
    res
  }
  # Return a tree with the gravity changed by factor
  .gravitise <- function (node) {
    # increase tree gravity
    # for root node, use root edge
    if (node == getSize (tree) + 1) {
      # get the following edges
      following.edges <- tree$edge[ ,1] == node
      # calculate amount to change the preceding and following edges
      change.by <- tree$root.edge * factor
      # add and remove by change.by
      tree$root.edge <- calc(tree$root.edge, change.by)
      for(following.edge in which(following.edges)) {
        new_following_edge <- calc(tree$edge.length[following.edge],
                                   change.by, FALSE)
        #if change has created negative branch length, add it to preceding edge length
        if(new_following_edge > 0) {
          tree$edge.length[following.edge] <- new_following_edge
        }
      }
      # add and remove by change.by
      tree$root.edge <- tree$root.edge + change.by
      tree$edge.length[following.edges] <-
        tree$edge.length[following.edges] - change.by
    } else {
      # get the following edges for node
      following.edges <- tree$edge[ ,1] == node
      if (any (following.edges)) {
        # get the preceding edges
        preceding.edge <- tree$edge[ ,2] == node
        # calculate amount to change the preceding and following edges
        change.by <- mean(tree$edge.length[following.edges]) * factor
        # add and remove by change.by
        tree$edge.length[preceding.edge] <- calc(tree$edge.length[preceding.edge],
                                                 change.by)
        for(following.edge in which(following.edges)) {
          new_following_edge <- calc(tree$edge.length[following.edge],
                                     change.by, FALSE)
          #if change has created negative branch length, add it to preceding edge length
          if(new_following_edge > 0) {
            tree$edge.length[following.edge] <- new_following_edge
          }
        }
      }
    }
    tree <<- tree
  }
  # make sure it is rooted
  if (!is.rooted (tree)) {
    stop ('Tree must be rooted!')
  }
  # make sure tree has root edge
  if (is.null (tree$root.edge)) {
    tree$root.edge <- 0
  }
  # list internal nodes
  internal.nodes <- length (tree$tip.label)+1:
    (length (tree$tip.label) +tree$Nnode)
  # if factor is negative, downGravitise
  if (factor < 0) {
    factor <- abs (factor)
    calc <- calcDwn
  } else {
    calc <- calcUp
  }
  plyr::m_ply (.data = data.frame (node = internal.nodes),
               .fun = .gravitise)
  tree
}

unbalancedTree <- function() {
  tree <- reBalance(rtree(n), -3)  # MTT function
  compute.brlen(tree)
}
balancedTree <- function() {
  tree <- reBalance(rtree(n), 3)
  compute.brlen(tree)
}
plotTreeWEds <- function() {
  ed_vals <- calcED(tree)
  tree$tip.label <- paste('t00000', 1:length(tree$tip.label))
  plot(tree, show.tip.label=TRUE, tip.color = "white")
  tiplabels(signif(ed_vals[,1], digits=2), adj = -.25, frame="none")
}

library(MoreTreeTools)
# plot parameters
n <- 12

tiff("~/Desktop/figure.tiff", width=9, height=9, units="cm",
     res=1200)
par(mfrow = c (2, 3), mar=c(1, 1, 1, 1))

# PF (balanced, slightly high gravity)
tree <- balancedTree()
for(i in 1:3) {
  # run in for loops with small factor to prevent negative branch lengths
  tree <- reGravitise(tree, .1)
}
plotTreeWEds()
# EMPTY
plot.new()
# EPH (balanced, slightly low gravity)
tree <- balancedTree()
for(i in 1:3) {
  tree <- reGravitise(tree, -.1)
}
plotTreeWEds()
# PAN (unbalanced, v. high gravity)
tree <- unbalancedTree()
for(i in 1:15) {
  tree <- reGravitise(tree, .1)
}
plotTreeWEds()
# HYD (unbalanced, mid high gravity)
tree <- unbalancedTree()
for(i in 1:10) {
  tree <- reGravitise(tree, .1)
}
plotTreeWEds()
# DE (balanced, high gravity)
tree <- unbalancedTree()
for(i in 1:3) {
  tree <- reGravitise(tree, .1)
}
plotTreeWEds()
dev.off()