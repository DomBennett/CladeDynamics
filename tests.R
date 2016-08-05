## 18/07/2014
## D.J. Bennett
## Confirm functions work

## Libs
library (testthat)
library (geiger)
source (file.path ('tools', 'test_tools.R'))
source (file.path ('tools', 'download_tools.R'))
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'parse_tools.R'))
source (file.path ('tools', 'precalculate_tools.R'))
source (file.path ('tools', 'simulation_tools.R'))

## Tests
cat ('\n\nRunning tests ...\n')
context ('Testing download tools ...')
test_that ('safeConnect([basic]) works ...', {
  expect_that (safeConnect (expr = a - 2,
                            wait = 0.1,
                            trys = 1,
                            verbose = FALSE),
               throws_error ())
  a <- 2
  expect_that (safeConnect (expr = a - 2,
                            wait = 0.1,
                            trys = 1,
                            verbose = FALSE),
               equals (0))
})
context ('Testing compare tools ...')
test_that ('calcTreeStats([basic]) works ...', {
  balanced.tree <- compute.brlen (stree (64, 'balanced'))
  unbalanced.tree <- compute.brlen (stree (64, 'left'))
  balanced.res <- calcTreeStats (list (balanced.tree))
  unbalanced.res <- calcTreeStats (list (unbalanced.tree))
  # colless should be lower for balanced tree
  expect_that (balanced.res[['colless']],
                    is_less_than(unbalanced.res[['colless']]))
})
context ('Testing parse tools ...')
test_that ('convertToDist([basic]) works ...', {
  # create a simple tree
  tree <- stree (10)
  res <- convertToDist (tree, 1)[[1]]
  expect_that (res$Nnode, equals (9))
})
test_that ('safeChronoMPL([basic]) works ...', {
  # create a simple non-ultrametric tree
  tree <- rtree (10)
  res1 <- safeChronoMPL (tree)
  expect_that (is.ultrametric (res1), is_true ())
  # create a simple tree without branch lengths
  tree <- stree (10)
  res2 <- safeChronoMPL (tree)
  expect_that (is.ultrametric (res2), throws_error ())
})
context ('Testing precalculate tools ...')
test_that ('pack([basic]) works ...', {
  # create a random mix of trees
  # ... too big
  res1 <- pack (stree (12, 'left'), 11, 9)
  expect_that (getSize (res1[[1]]), equals (11))
  # ... too small
  res2 <- pack (stree (8, 'left'), 11, 9)
  expect_that (res2, is_null ())
  # ... some big, some small
  trees <- list ()
  ns <- c (8,8,8,8,8,12,12,12,12,12)
  for (i in 1:length (ns)) {
    n <- ns[i]
    trees <- c (trees, list (stree (n, 'left')))
  }
  class (trees) <- 'multiPhylo'
  res3 <- pack (trees, 11, 9)
  expect_that (length (res3), equals (5))
})
context ('Testing EDBMM model')
test_that ('.seedTree([basic]) works ...', {
  tree <- .seedTree (10, 10)
  # should have ten tips
  expect_that (length (tree$tip.label), equals (10))
  # should be 10 time units old
  expect_that (getAge (tree, node=length (tree$tip.label) + 1)$age,
               equals (10))
})
test_that ('runEDBMM([basic]) works ...', {
  # grow a tree by 10 tips
  tree <- runEDBMM (birth=2, death=1, stop.at=10, fossils=TRUE)
  # it should have 12 tips (because a default seed of two starts it)
  expect_that (getSize (drop.extinct(tree)), equals (12))
  # it should have extinct tips
  cat (paste0 ('NOTE! There is a 0.017 probability of \'runEDBMM([basic])\' test failing.
               Try running again if it does fail.'))
  expect_that (is.ultrametric (tree), is_false ())
})
test_that ('runEDBMM(record=TRUE) works ...', {
  # grow a tree by 20 tips, record tree growth every 10 tips added
  tree <- runEDBMM (birth=1, death=0.5, stop.at=20, sample.at=10,
                    record=TRUE)
  expect_that (getSize (tree[[1]]), equals (2))
  expect_that (getSize (tree[[2]]), equals (12))
  expect_that (getSize (tree[[3]]), equals (22))
})
cat ('\n\nTests passed.\n')