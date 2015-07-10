## 29/08/2014
## D.J. Bennett
## Unused functions created as part of pipeline

plotWithNodeLabels <- function (tree) {
  plot (tree, no.margin = TRUE)
  nodelabels (text = tree$node.label)
}

## Creates a nice little gif of a tree through time, needs imagemgik
plotTreeGrowth <- function (trees, file.dir, time.steps = FALSE,
                            animation.time = 10){
  ## Generate a .gif showing a tree grow
  if (time.steps[1] == FALSE) {
    time.steps <- 1:length (trees)
  }
  delay <- animation.time/length (trees)
  par (mar = c (2, 2, 3, 2) + 0.1)
  png (file = "temp_tree_%010d.png")
  last.tree <- trees[[length (trees)]]
  x.lim <- getAge (last.tree, length (last.tree$tip.label) + 1)
  for (i in 1:length (trees)) {
    plot (trees[[i]], show.tip.label = FALSE, x.lim = x.lim)
    mtext (paste0 ('t = ', time.steps[i]), adj = 0.1, line = 0)
  }
  closeDevices ()
  system (paste0 ("convert -delay ", delay," *.png ", file.dir))
  file.remove (list.files (pattern = ".png"))
}

## For plotting fates against ED, would need to keep extinct species in tree
##  this would make the pipeline take longer to run -- not worth it
getFates <- function (trees, sample) {
  ## Determine the fate of every species in all time steps
  spp <- paste0 ('t', 1:length (trees[[length (trees)]]$tip.label))
  # create a new list of trees with extinct species dropped
  extant.trees <- list (trees[[1]]) # first tree is the seed, everything is extant
  .dropExtinct <- function (i) {
    tree <- drop.fossil (trees[[i]])
    extant.trees <<- c (extant.trees, list (tree))
  }
  m_ply (.data = data.frame (i = 2:length (trees)), .fun = .dropExtinct)
  # workout which species were extant and went extinct using new list
  .eachTree <- function (i) {
    last.tree <- extant.trees[[i-1]]
    this.tree <- extant.trees[[i]]
    extant <- this.tree$tip.label
    extinct <- last.tree$tip.label[!last.tree$tip.label %in%
                                     this.tree$tip.label]
    survived <- last.tree$tip.label[last.tree$tip.label %in%
                                      this.tree$tip.label]
    children <- this.tree$tip.label[!this.tree$tip.label %in%
                                      last.tree$tip.label]
    # parents are those that survived and have a tip.edge equal or less than
    #  children lengths
    parents <- c ()
    max.child.length <- max (this.tree$edge.length[
      this.tree$edge[ ,2] %in% which (this.tree$tip.label %in% children)])
    for (each in survived) {
      each.length <- this.tree$edge.length[
        this.tree$edge[ ,2] %in% which (this.tree$tip.label == each)]
      if (each.length <= max.child.length) {
        parents <- c (parents, each)
      }
    }
    list (extinct = extinct, extant = extant, survivors = survived,
          speciators = parents)
  }
  species.success <-
    mlply (.data = data.frame (i = 2:length (trees)), .fun = .eachTree)
  # create data frame of success
  res <- matrix (ncol = length (species.success), nrow = length (spp))
  rownames (res) <- spp
  .eachRow <- function (i) {
    res[species.success[[i]]$extinct, i] <<- -1
    res[species.success[[i]]$survivors, i] <<- 0
    res[species.success[[i]]$speciators, i] <<- 1
  }
  m_ply (.data = data.frame (i = 1:length (species.success)),
         .fun = .eachRow)
  res
}

getEDs <- function (trees) {
  eds <- list ()
  .calc <- function (i) {
    res <- calcED (trees[[i]])
    eds <<- c (eds, list (res))
  }
  m_ply (.data = data.frame (i = 1:length (trees)),
         .fun = .calc)
  eds
}

plotFateVsED <- function (fates, eds, time.lag = 1) {
  # find corresponding x and y for each time + lag
  x <- y <- c ()
  for (i in (time.lag + 1):(length (trees) - time.lag)) {
    fate.slice <- fates[!is.na (fates[ ,i]),i]
    ed.slice <- eds[[i-time.lag]]
    x <-
      c (x, na.omit (ed.slice[match (names (fate.slice), ed.slice[ ,1]),2]))
    y <-
      c (y, na.omit (fate.slice[match (ed.slice[ ,1], names (fate.slice))]))
  }
  # plot
  plot (x = x, y = y, xlab = 'ED', ylab = 'Species\' Fate',
        col = rainbow (3, alpha = 0.5)[3], pch = 19,
        main = paste0 ('Time lag: [', time.lag, ']'))
}

# test_that ('getFates([basic]) works ...', {
#   trees <- genRandomTrees (100)
#   # fates should all be 0 and 1 at the beginning
#   # because there is no extinction in the random trees
#   # and species are added to the tree at the newset tip
#   fates <- getFates (trees)
#   # so... the tree starts with 3 tips, so the first fate for 
#   # tip 3 is 1, this is where tip 4 is added in time step 2
#   expect_that (as.numeric (fates [3,1]), equals (1))
# })
# test_that ('plotFatesVsED([basic]) works ...', {
#   trees <- genRandomTrees (10)
#   fates <- getFates (trees)
#   eds <- getEDs (trees)
#   plotFateVsED (fates, eds, time.lag = 1)
#   expect_that (is.null (dev.list ()), is_false ())
#   dev.off ()
# })

## Defunct functions for identifying the rise and fall of a clade, for plotting
##  Now I might use, a spindle diagram? Or use CM to identify the centre?
# findRiseAndFall <- function (element, min.size, min.time) {
#   ## Find a radiation by identifying the peak and local troughs
#   ## Be selective by choosing only those clades that have a peak
#   ## above min.size, and read their troughs above min.time
#   # only extinct clades
#   if (element[length (element)] != 0) {
#     return (NA)
#   }
#   if (max (element) < min.size) {
#     return (NA)
#   }
#   # find peak (only use first)
#   peak <- which (element == max (element))[1]
#   # find elements that increase up to peak and decrease after, set peak as false
#   bool <- c (element[1:(peak - 1)] <= element[2:peak], FALSE,
#      element[peak:(length (element) - 1)] >= element[(peak + 1):length (element)])
#   # identify start and end by finding immediate falses around peak
#   falses <- which (bool == FALSE)
#   if (falses[1] == peak) {
#     # if first of the falses is peak, the first element is the start
#     start <- 1
#   } else {
#     start <- falses[which (falses == peak) - 1] + 1
#   }
#   end <- falses[which (falses == peak) + 1] - 1
#   if (length (start:end) < min.time) {
#     return (NA)
#   }
#   radiation.pos <- start:end
# }

findRiseAndFall <- function (element, min.size, min.time) {
  ## Find peak and centre time around it
  if (element[length (element)] != 0) {
    return (NA)
  }
  if (max (element) < min.size) {
    return (NA)
  }
  # find peak (only use first)
  peak <- which (element == max (element))
  
  # find elements that increase up to peak and decrease after, set peak as false
  bool <- c (element[1:(peak - 1)] <= element[2:peak], FALSE,
             element[peak:(length (element) - 1)] >= element[(peak + 1):length (element)])
  # identify start and end by finding immediate falses around peak
  falses <- which (bool == FALSE)
  if (falses[1] == peak) {
    # if first of the falses is peak, the first element is the start
    start <- 1
  } else {
    start <- falses[which (falses == peak) - 1] + 1
  }
  end <- falses[which (falses == peak) + 1] - 1
  if (length (start:end) < min.time) {
    return (NA)
  }
  radiation.pos <- start:end
}

# findRiseAndFall <- function (element, min.size = 2) {
#   ## Find all clade radiations that rise from and to min.size
#   bool <- element >= min.size
#   # add shifted bools to get a 1,2 or 3 representing the
#   #  position of the number in the sequence
#   n.pos <- bool + c (bool[-1], bool[1]) +
#     c (bool[length (bool)], bool[-length (bool)])
#   # find pos 1s that mark the beginning and end of the radiation
#   # with regexp
#   n.pos.str <- paste (n.pos, collapse = '')
#   start <- regexpr ('123', n.pos.str)
#   end <- regexpr ('321', n.pos.str)
#   if (any (c (start, end) < 0)) {
#     return (NA)
#   }
#   # return radiation index
#   start:(end+2)
# }

## A way to plot range of ages of clades
histClageAges <- function (metadata) {
  ## Take metadata and plot clade ages in hist for each ED strength
  ages <- strength <- NULL
  for (i in 1:nrow (metadata)) {
    res <- read.csv (file = file.path (
      'results', metadata$res.dir[i], 'clades_through_time.csv'))[ ,-1]
    ages <- c (ages, colSums (res != 0))
    strength <- c (strength, rep (metadata$strength[i], ncol (res)))
  }
  clade.ages <- data.frame (ages = ages, strength = strength)
  hist <- ggplot (clade.ages, aes (x = ages, fill = strength))
  hist + geom_density ()
}

# if record, plot clade stats against ED strength
if (record) {
  ed.strengths <- timespans <- cgs <- cms <- NULL
  cladestatsfiles <- sub ('\\.tre', '\\.csv',
                          metadata$treefilename)
  for (i in 1:nrow (metadata)) {
    clade.stats <- read.csv (file.path (res.dir, cladestatsfiles[i]))
    timespans <- c (timespans, clade.stats$time.span)
    cgs <- c (cgs, clade.stats$cg)
    cms <- c (cms, clade.stats$cm)
    ed.strengths <- c (ed.strengths,
                       rep (metadata$strength[i], nrow (clade.stats)))
  }
  pdf (file.path (res.dir, 'cladestats_ED_strength.pdf'))
  plot (timespans ~ ed.strengths, xlab = 'ED strength', ylab = 'Clade time span',
        col = rainbow (3, alpha = 0.7)[3], pch = 19)
  plot (cms ~ ed.strengths, xlab = 'ED strength', ylab = 'Centre of Mass',
        col = rainbow (3, alpha = 0.7)[3], pch = 19)
  plot (cgs ~ ed.strengths, xlab = 'ED strength', ylab = 'Centre of gyration',
        col = rainbow (3, alpha = 0.7)[3], pch = 19)
  closeDevices ()
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
  pdf (file.path (res.dir, 'pca.pdf'))
  # remove any that aren't ultrametric or rate.smooted
  real.stats <- real.stats[real.stats$ultra | real.stats$chronos, ]
  real.stats <- real.stats[real.stats$ultra, ]
  real.stats$sig <- NA
  real.stats$eps <- NA
  cols <- c ('sig', 'eps', 'colless', 'sackin', 'tci', 'gamma')
  input <- rbind (stats[ ,cols], real.stats[, cols])
  pca.res <- prcomp (input[,cols[-c(1,2)]],
                     scale. = TRUE, center = TRUE)
  pca.x <- as.data.frame(pca.res$x[!is.na (input$eps), ])
  pca.x.real <- as.data.frame(pca.res$x[is.na (input$eps), ])
  # get pca distances from real
  quick <- function (x, real) {
    min (abs (x - real))
  }
  # use min dist
  stats$pc1.dists <- sapply (pca.x$PC1, quick, pca.x.real$PC1)
  stats$pc2.dists <- sapply (pca.x$PC2, quick, pca.x.real$PC2)
  stats$pc3.dists <- sapply (pca.x$PC3, quick, pca.x.real$PC3)
  # use mean dist
  stats$pc1.dists <- abs (pca.x$PC1 - mean.real['PC1'])
  stats$pc2.dists <- abs (pca.x$PC2 - mean.real['PC2'])
  stats$pc3.dists <- abs (pca.x$PC3 - mean.real['PC3'])
  stats$PC1 <- pca.x$PC1
  stats$PC2 <- pca.x$PC2
  stats$PC3 <- pca.x$PC3
  pca.rot <- as.data.frame (pca.res$rotation)
  prop.var <- round (sapply (pca.res$sdev^2,
                             function (x) Reduce('+', x)/sum (pca.res$sdev^2)), 3)
  names (prop.var) <- colnames (pca.rot)
  comparisons <- list (c ("PC1", "PC2"), c ("PC2", "PC3"), c ("PC1", "PC3"))
  for (comp in comparisons) {
    psize <- 10
    p <- ggplot (pca.x, aes_string (x = comp[1], y = comp[2])) +
      geom_point (aes (colour = 'blue'), size = psize) +
      #scale_colour_gradient2 (low = "red", high = "blue", name = expression(psi),
      #                        guide = guide_legend(keywidth = 5, keyheight = 5)) +
      geom_point (data = pca.x.real, colour = 'black', shape = 15, size = psize) +
      xlab (paste0 (comp[1], " (", prop.var[comp[1]], ")")) +
      ylab (paste0 (comp[2], " (", prop.var[comp[2]], ")")) +
      theme_bw (base_size=48)
    print (p)
    rm (p)
    plot(x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]], xlab = comp[1],
         ylab = comp[2], cex = 0.5, pch = 19)
    text (x = pca.rot[ ,comp[1]], y =  pca.rot[ ,comp[2]],
          rownames (pca.rot[comp[1]]), adj = 1)
  }
  closeDevices ()
  return (pca.dists)
}

# context ('Testing analysis tools')
# test_that ('plotSuccess([basic]) works ...', {
#   # generate random clade data using Poisson dist
#   data <- genRandomData (100, 100, 1)
#   # plot
#   plotSuccess (data)
#   # if plot is generated check devices
#   expect_that (is.null (dev.list ()), is_false ())
#   # turn off devices
#   dev.off ()
# })
# test_that ('plotNormalisedSuccess([basic]) works ...', {
#   # generate random clade data using Poisson dist
#   data <- genRandomData (100, 100, 1)
#   # plot
#   plotNormalisedSuccess (data)
#   # if plot is generated check devices
#   expect_that (is.null (dev.list ()), is_false ())
#   # turn off devices
#   dev.off ()
# })
# test_that ('.reformat([basic]) works ...', {
#   # generate random clade success data using Poisson dist.
#   clade.success <- list ()
#   for (i in 1:10) {
#     node <- paste0 ('n', unique (rpois (100, i)))
#     n.children <- rpois (length (node), 10 + i)
#     clade.success <- c (clade.success, list (
#       data.frame (node = node, n.children = n.children)))
#   }
#   reformatted <- .reformat (clade.success, 1)
#   # test that one of the elements is as expected
#   node <- as.character (clade.success[[1]][1,1])
#   n.children <- clade.success[[1]][1,2]
#   expect_that (reformatted[node][1,], equals (n.children))
# })
# test_that ('.countChildren([basic]) works ...', {
#   # simple tree of 10 tips and one internal node
#   tree <- stree (10)
#   # .countChildren requires node.labels
#   tree$node.label <- 'n1'
#   # .countChildren requires a list of extinct species
#   extinct <- c ()
#   res1 <- .countChildren (tree, extinct)
#   expect_that (res1[1,2], equals (10))
#   # adding a species to extinct ...
#   extinct <- c ('t1')
#   res2 <- .countChildren (tree, extinct)
#   expect_that (res2[1,2], equals (9))
# })
# ## TODO:
# ##  --calcCladeSuccess
# ##  --getCladeSuccess