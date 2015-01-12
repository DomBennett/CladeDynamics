
library (MoreTreeTools)
library (apTreeshape)


# explore fractal properties of EDBMM model

# r: births and deaths
# n: starting tree size


shapeModel <- function (r, n, t) {
  addTips <- function(target) {
    # add tip
    new.tip <- paste0 ('t', tip.counter)
    edge <- which (tree$edge[ ,2] == which (tree$tip.label == target))
    tree <- addTip (tree = tree, edge = edge,
                    tip.name = new.tip, node.age = 0)
    tip.counter <<- tip.counter + 1
    tree
  }
  run <- function (i) {
    # calc ED
    eds <- calcED (tree)
    eds.order <- order (eds[,1])
    # identify growth and shrink spots
    growth.spots <- rownames (eds)[eds.order[1:r]]
    shrink.spots <- rownames (eds)[eds.order[(nrow (eds)-(r-1)):nrow (eds)]]
    # add and remove accordingly
    tree <<- drop.tip (tree, tip = shrink.spots)
    for (target in growth.spots) {
      tree <<- addTips (target)
    }
    # grow
    edges <- which (tree$edge[ ,2] <= length (tree$tip.label))
    tree$edge.length[edges] <- tree$edge.length[edges] + 1/n
    tree <<- tree
  }
  if (r > n/2) {
    stop ('Error: max r is n/2')
  }
  tip.counter <- n + 1
  tree <- compute.brlen (rtree (n))
  tree$edge.length <- tree$edge.length/getSize (tree)
  m_ply (.data = data.frame (i = 1:t), .fun = run)
  tree
}

bifurcationModel <- function (rs, ns, ts, iterations=100) {
  runmodel <- function (r, n, t) {
    cat ('.... n = [', n, ']', sep = '')
    tree <- shapeModel (r, n, t)
    ts.tree <- as.treeshape (tree)
    colless (ts.tree, 'yule')
  }
  rs <- rep (rs, each = iterations)
  mdply (.data = data.frame (r=rs, n=ns, t=ts), .fun = runmodel)
}

rs <- 1:6
ns <- rep (100, length (rs))
ts <- rep (100, length (rs))
res <- bifurcationModel(rs=rs, ns=ns, ts=ts, iterations = 20)
plot (res$V1~res$r, ylab = 'Imbalance', xlab = 'r')
res

# NOT CHAOTIC