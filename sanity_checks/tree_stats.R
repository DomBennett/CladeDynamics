## 24/07/2014
## D.J. Bennett
## Check if my measures for controlling different tree sizes works

## Libs
source (file.path ('tools', 'compare_tools.R'))

## Create different sizes of trees using the same process -- birth-death
how.many.trees <- 1000
ns <- runif (how.many.trees, min = 100, max = 1000)
trees <- list ()
modelTrees <- function (i) {
  tree <- drop.extinct (sim.bdtree (b = 2, d = 1, n = ns[i],
                                    stop = 'taxa', extinct = FALSE))
  trees <<- c (trees, list (tree))
}
m_ply (.data = data.frame (i = 1:how.many.trees),
       .fun = modelTrees, .progress = 'time')

## Test their stats ...
res <- calcMeanTreeShapeStats (trees)
# all should be normally distributed
hist (res$colless.stat)
hist (res$sackin.stat)
hist (res$iprime.stat)
hist (res$tc.stat)
hist (res$gamma.stat)
# there should be no relationship between size and stat
plot (x = ns, y = res$colless.stat, col = 'cornflowerblue', pch = 19)
abline (lm (res$colless.stat ~ ns), col = 'red')
plot (x = ns, y = res$sackin.stat, col = 'cornflowerblue', pch = 19)
abline (lm (res$sackin.stat ~ ns), col = 'red')
plot (x = ns, y = res$iprime.stat, col = 'cornflowerblue', pch = 19)
abline (lm (res$iprime.stat ~ ns), col = 'red')
plot (x = ns, y = res$tc.stat, col = 'cornflowerblue', pch = 19)
abline (lm (res$tc.stat ~ ns), col = 'red')
plot (x = ns, y = res$gamma.stat, col = 'cornflowerblue', pch = 19)
abline (lm (res$gamma.stat ~ ns), col = 'red')
# All looks as expected. Phew.... sanity reconfirmed!