## 18/07/2014
## D.J. Bennett
## Confirm functions work

## Libs
library (testthat)
source (file.path ('tools', 'test_tools.R'))
source (file.path ('tools', 'model_tools.R'))
source (file.path ('tools', 'analysis_tools.R'))
source (file.path ('tools', 'compare_tools.R'))

## Tests
# Model tools first ...
context ('Testing model tools ...')
test_that ('seedTree([basic]) works ...', {
  tree <- seedTree (10, 10)
  # should have ten tips
  expect_that (length (tree$tip.label), equals (10))
  # should be 10 time units old
  expect_that (getAge (tree, length (tree$tip.label) + 1),
               equals (10))
})
test_that ('.countChildren([basic]) works ...', {
  # simple tree of 10 tips and one internal node
  tree <- stree (10)
  # .countChildren requires node.labels
  tree$node.label <- 'n1'
  # .countChildren requires a list of extinct species
  extinct <- c ()
  res1 <- .countChildren (tree, extinct)
  expect_that (res1[1,2], equals (10))
  # adding a species to extinct ...
  extinct <- c ('t1')
  res2 <- .countChildren (tree, extinct)
  expect_that (res2[1,2], equals (9))
})
test_that ('.reformat([basic]) works ...', {
  # generate random clade success data using Poisson dist.
  clade.success <- list ()
  for (i in 1:10) {
    node <- paste0 ('n', unique (rpois (100, i)))
    n.children <- rpois (length (node), 10 + i)
    clade.success <- c (clade.success, list (
      data.frame (node = node, n.children = n.children)))
  }
  reformatted <- .reformat (clade.success, 1)
  # test that one of the elements is as expected
  node <- as.character (clade.success[[1]][1,1])
  n.children <- clade.success[[1]][1,2]
  expect_that (reformatted[node][1,], equals (n.children))
})
test_that ('growMRMMTree([basic]) works ...', {
  # grow a tree to 10 tips
  tree <- growMRMMTree (birth = 2, death = 1, stop.at = 10)
  # it should have 10 tips
  expect_that (length (drop.extinct(tree)$tip.label), equals (10))
  # it should have extinct tips
  cat (paste0 ('NOTE! There is a 0.017 probability of \'growMRMMTree([basic])\' test failing.
               Try running again if it does fail.'))
  expect_that (is.ultrametric (tree), is_false ())
})
# Now for those analysis tools ...
context ('Testing analysis tools')
test_that ('plotSuccess([basic]) works ...', {
  # generate random clade data using Poisson dist
  data <- genRandomData (100, 100, 1)
  # plot
  plotSuccess (data)
  # if plot is generated check devices
  expect_that (is.null (dev.list ()), is_false ())
  # turn off devices
  dev.off ()
})
test_that ('plotNormalisedSuccess([basic]) works ...', {
  # generate random clade data using Poisson dist
  data <- genRandomData (100, 100, 1)
  # plot
  plotNormalisedSuccess (data)
  # if plot is generated check devices
  expect_that (is.null (dev.list ()), is_false ())
  # turn off devices
  dev.off ()
})
test_that ('getFates([basic]) works ...', {
  trees <- genRandomTrees (100)
  # fates should all be 0 and 1 at the beginning
  # because there is no extinction in the random trees
  # and species are added to the tree at the newset tip
  fates <- getFates (trees)
  # so... the tree starts with 3 tips, so the first fate for 
  # tip 3 is 1, this is where tip 4 is added in time step 2
  expect_that (as.numeric (fates [3,1]), equals (1))
})
test_that ('plotFatesVsED([basic]) works ...', {
  trees <- genRandomTrees (10)
  fates <- getFates (trees)
  eds <- getEDs (trees)
  plotFateVsED (fates, eds, time.lag = 1)
  expect_that (is.null (dev.list ()), is_false ())
  dev.off ()
})
# untested (not critical to model):
#  nomarlise (sub function of plotNormalisedSuccess)
#  findRiseAndFall (not used anymore)
#  findNonZeros (sub function of plotSuccess)
#  plotTreeGrowth (using an external program + subject to change)
context ('Testing compare tools ...')
test_that ('.calcTreeShapeStats ([basic]) works ...', {
  test.tree <- stree (64, 'balanced')
  test.tree <- compute.brlen (test.tree)
  res <- .calcTreeShapeStats (test.tree)
  # colless test should be 0 for a balanced tree
  expect_that (res[['colless.stat']], equals (0))
})
test_that ('calcTreeShapeStates([basic]) works ...', {
  test.trees <- list ()
  for (i in 1:10) {
    test.tree <- stree (64, 'balanced')
    test.tree <- compute.brlen (test.tree)
    test.trees <- c (test.trees, list (test.tree))
  }
  res <- calcTreeShapeStats (test.trees)
  # colless test should be 0 for balanced trees
  expect_that (res[['mean.colless.stat']], equals (0))
})
test_that ('extractStat([basic]) works ...', {
  # create test simulated list
  simulated.tree.stats <-
    list ('1' = list (mean = 1, sd = 1, stat = c (1,1,1,1)),
          '0' = list (mean = 1, sd = 1, stat = c (1,1,1,1)),
          '-1' = list (mean = 1, sd = 1, stat = c (1,1,1,1)))
  res1 <- extractStat (simulated.tree.stats, 'mean')
  expect_that (res1, equals (c (1,1,1)))
  res2 <- extractStat (simulated.tree.stats, 'stat')
  expect_that (unlist (res2), equals (rep (1, 12)))
})
# untested (not critical to model):
#  histCladeAge (requires writing test files, likely to change)
cat ('\n\nTests passed.')
