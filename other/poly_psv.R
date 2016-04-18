# Do npolytomies correlate with PSV?

source (file.path ('tools', 'compare_tools.R'))
data.dir <- file.path ('data', 'treestats')
empirical.filenames <- c ('min50_max500_rspathD8.Rd',
                          'min50_max500_rschronopl.Rd',
                          'min50_max500_rschronoMPL.Rd')
# load pre-calculated empirical tree stats -- real.stats
real.stats <- readInMultiple (file.path (data.dir, empirical.filenames))

# load tree and count polys
pply <- rep(NA, nrow(real.stats))
for (i in 1:nrow(real.stats)) {
  print(i)
  fn <- as.character(real.stats$filename[i])
  tree <- read.tree(file.path('data', 'parsed_trees', fn))
  if(is(tree, 'phylo')) {
    pply[i] <- sum(tree$edge.length == 0)/length(tree$edge.length)
  } else {
    tpply <- 0
    for(t in tree) {
      tpply <- tpply + sum(t$edge.length == 0)/length(t$edge.length)
    }
    pply[i] <- tpply/length(tree)
  }
}

psv <- real.stats$psv.pathD8[!is.na(real.stats$psv.pathD8)]
pply2 <- pply[!is.na(real.stats$psv.pathD8)]
plot(psv~pply2)
cor.test(psv, pply2)  # the more polytomies the higher the gravity

tapply(pply, real.stats$quality, mean)
tapply(real.stats$psv.pathD8, real.stats$quality, mean)
