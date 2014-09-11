## 18/07/2014
## D.J. Bennett
## Confirm functions work

## Libs
library (testthat)
source (file.path ('tools', 'test_tools.R'))
source (file.path ('tools', 'model_tools.R'))
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'parse_tools.R'))
source (file.path ('tools', 'precalculate_tools.R'))

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
context ('Testing compare tools ...')
test_that ('calcTreeShapeStates([basic]) works ...', {
  test.trees <- list ()
  for (i in 1:5) {
    test.tree <- stree (64, 'balanced')
    test.tree <- compute.brlen (test.tree)
    test.trees <- c (test.trees, list (test.tree))
  }
  res <- calcTreeShapeStats (test.trees, iterations = 1)
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
context ('Testing parse tools ...')
test_that ('convertToDist([basic]) works ...', {
  # create a simple tree
  tree <- stree (10)
  res <- convertToDist (tree, 1)[[1]]
  expect_that (res$Nnode, equals (9))
})
test_that ('safeChronos([basic]) works ...', {
  # create a simple non-ultrametric tree
  tree <- rtree (10)
  res1 <- safeChronos (tree, lambda = 1, quiet = TRUE)
  expect_that (is.ultrametric (res1), is_true ())
  # create a simple tree without branch lengths
  tree <- stree (10)
  res2 <- safeChronos (tree, lambda = 1, quiet = TRUE)
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
cat ('\n\nTests passed.')