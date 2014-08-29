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
test_that ('runEDBMM([basic]) works ...', {
  # grow a tree by 10 tips
  tree <- runEDBMM (birth = 2, death = 1, stop.at = 10, fossils = TRUE)
  # it should have 12 tips (because a default seed of two starts it)
  expect_that (getSize (drop.extinct(tree)), equals (12))
  # it should have extinct tips
  cat (paste0 ('NOTE! There is a 0.017 probability of \'runEDBMM([basic])\' test failing.
               Try running again if it does fail.'))
  expect_that (is.ultrametric (tree), is_false ())
})
test_that ('runEDBMM(record=TRUE) works ...', {
  # grow a tree by 20 tips, record tree growth every 10 tips added
  tree <- runEDBMM (birth = 2, death = 1, stop.at = 20, sample.at = 10,
                    record = TRUE)
  expect_that (getSize (tree[[1]]), equals (2))
  expect_that (getSize (tree[[2]]), equals (12))
  expect_that (getSize (tree[[3]]), equals (22))
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
## TODO:
##  --calcCladeSuccess
##  --getCladeSuccess
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
cat ('\n\nTests passed.')