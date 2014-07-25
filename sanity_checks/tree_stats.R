## 24/07/2014
## D.J. Bennett
## Check if my measures for controlling different tree sizes works

## Libs
source (file.path ('tools', 'compare_tools.R'))

## Create different sizes of trees using the same process -- ERMM
how.many.trees <- 100
min.ntips <- 10
stop.by <- sample (c ('taxa', 'time'), how.many.trees, replace = TRUE)
ns <- runif (how.many.trees, min = 100, max = 5000)
trees <- list ()
for (i in 1:length (ns)) {
  tree <- sim.bdtree (b = 2, d = 1, n = ns[i],
                      stop = stop.by[i], extinct = FALSE)
  if (length (tree$tip.label) > min.ntips) {
    # avoid using really small trees
    trees <- c (trees, list (tree))
  }
}

## Test their stats ...
res <- calcMeanTreeShapeStats (trees)
# all should be normally distributed
hist (res$colless.stat)
hist (res$sackin.stat)
hist (res$iprime.stat)
hist (res$tc.stat)
hist (res$gamma.stat)