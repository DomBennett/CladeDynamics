## 16/07/2014
## D.J. Bennett
## Comparing different model runs and natural trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))

## Parameters
if (!exists ('pars')) {
  pars <- list (n.model = 2, seed = 2,
                max.birth = 5, min.birth = 1.1,
                max.death = 1, min.death = 1,
                bias = 'FP', stop.by = 'n',
                max.ntaxa = 200, min.ntaxa = 50,
                min.psi = -1, max.psi = 1,
                reference = TRUE, iterations = 100)
  name <- 'testset'
}

## Dirs
data.dir <- file.path ('data', 'treestats')

## Functions
readTrees <- function (metadata, res.dir, runlog) {
  trees <- list ()
  for (i in 1:nrow (metadata)) {
    # read in tree
    tree.dir <- file.path (res.dir, metadata$treefilename[i])
    tree <- read.tree (tree.dir)
    # make sure it isn't multiPhylo
    if (class (tree) == 'multiPhylo') {
      tree <- tree[[length (tree)]]
    }
    # then add to list
    trees <- c (trees, list (tree))
  }
  trees
}

calcStats <- function (trees, metadata, iterations, reference,
                       res.dir) {
  # loop through each tree
  tree.stats <- list ()
  for (i in 1:length (trees)) {
    temp.res <- list (calcTreeShapeStats (
      trees[i], iterations = iterations, reference = reference))
    tree.stats <- c (tree.stats, temp.res)
  }
  # write out
  save (tree.stats, file = file.path (res.dir, 'shapestats.Rd'))
  # TODO: rethink how the calcTreeShapeStats works to avoid this loop
  colless.stat <- sackin.stat <- iprime.stat <- 
    gamma.stat <- tc.stat <- rep (NA, length (tree.stats))
  for (i in 1:length (tree.stats)) {
    colless.stat[i] <- tree.stats[[i]][['colless.stat']]
    sackin.stat[i] <- tree.stats[[i]][['sackin.stat']]
    iprime.stat[i] <- tree.stats[[i]][['iprime.stat']]
    gamma.stat[i] <- tree.stats[[i]][['gamma.stat']]
    tc.stat[i] <- tree.stats[[i]][['tc.stat']]
  }
  stats <- data.frame (colless.stat, sackin.stat, iprime.stat,
                       gamma.stat, tc.stat)
  stats <- cbind (metadata, stats)
  stats
}

compare <- function (stats, real.stats) {
  res <- data.frame ()
  window.size <- 0.5
  stat.names <- c ('colless', 'sackin', 'gamma', 'tci')
  for (i in 1:length (stat.names)) {
    diff <- windowAnalysis (stats[ ,stat.names[i]],
                            real.stats[ ,stat.names[i]],
                            size = window.size,
                            psi = stats[, 'psi'], incr = 0.1)
    temp <- data.frame (diff, size = window.size,
                        stat = stat.names[i])
    res <- rbind (res, temp)
  }
  res <- na.omit (res)
  res
}

plotResults <- function (res, stats, res.dir, metadata) {
  pdf (file.path (res.dir, 'treestats_ED_strength.pdf'))
  # plot distributions differences
  p <- ggplot (res, aes (mid, diff))
  print (p + geom_point (aes (colour = stat)))
  # plot psi and stat
  stat.names <- c ('colless', 'sackin', 'gamma', 'tci')
  par (mfrow = c (2,2), mar = c (3, 5, 0.5, 0.5))
  y <- metadata$psi
  for (each in stat.names) {
    plot (x = stats[ ,each], y = y, ylab = expression (psi), pch = 19,
          col = rainbow (3, alpha = 0.8)[3], cex.lab = 3, cex = 3, xlab = '',
          cex.axis = 2)
    model <- lm (y ~ stats[ ,each])
    abline (model)
  }
  closeDevices ()
}

pca <- function (stats, real.stats) {
  input <- rbind (stats, real.stats)
  pca.res <- prcomp (input[ ,c ('gamma', 'tci', 'sackin', 'colless')],
                     scale. = TRUE, center = TRUE)
  pca.x <- as.data.frame(pca.res$x[!is.na (input$psi), ])
  pca.x.real <- as.data.frame(pca.res$x[is.na (input$psi), ])
  pca.rot <- as.data.frame (pca.res$rotation)
  prop.var <- round (sapply (pca.res$sdev^2,
                             function (x) Reduce('+', x)/sum (pca.res$sdev^2)), 3)
  names (prop.var) <- colnames (pca.rot)
  comparisons <- list (c ("PC1", "PC2"), c ("PC2", "PC3"), c ("PC1", "PC3"))
  for (comp in comparisons) {
    p <- ggplot (pca.x, aes_string (x = comp[1], y = comp[2])) +
      geom_point (aes (colour = input$psi[!is.na (input$psi)])) +
      scale_colour_gradient2 (low = "red", high = "blue", name = expression(psi)) +
      geom_point (data = pca.x.real, colour = 'black', shape = 15) +
      xlab (paste0 (comp[1], " (", prop.var[comp[1]], ")")) +
      ylab (paste0 (comp[2], " (", prop.var[comp[2]], ")")) +
      theme_bw ()
    print (p)
    rm (p)
    plot(x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]], xlab = comp[1],
         ylab = comp[2], cex = 0.5, pch = 19)
    text (x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]],
          rownames (pca.rot[comp[1]]), adj = 1)
  }
}


## Run
res.dir <- file.path ('results', name)
runlog <- file.path (res.dir, 'runlog.csv')
cat ('\nReading in data ...')
# get metadata
metadata <- read.csv (runlog, stringsAsFactors = FALSE)
# load pre-calculated natural tree stats
filename <- paste0 ('min', pars$min.ntaxa, '_max', pars$max.ntaxa, '.Rd')
load (file.path (data.dir, filename))
trees <- readTrees (metadata, res.dir, runlog)
cat ('\nCalculating tree stats ...')
stats <- calcTreeStats(trees)
stats <- cbind (metadata, stats)
cat ('\nCompare to real trees ...')
res <- compare (stats, real.stats)
cat ('\nPlotting ...')
plotResults (res, stats, res.dir, metadata)
pca (stats, real.stats)