## 18/07/2014
## D.J. Bennett
## Confirm functions work

## Libs
library (testthat)
source (file.path ('tools', 'model_tools.R'))
source (file.path ('tools', 'analysis_tools.R'))

## Tests
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
  expect_that (length (tree$tip.label), equals (10))
})