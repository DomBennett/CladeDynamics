## 17/07/2014
## D.J. Bennett
## Tools for comparing MRMM tree runs

## Libs
library (MoreTreeTools)
library (ape)
library (apTreeshape)
library (caper)
library (geiger)
library (ggplot2)

pca2 <- function (stats, real.stats, stat.names, filename,
                 ignore.chronos=TRUE) {
  addDistances <- function (pc.sim) {
    getDist <- function (exp.sig, exp.eps) {
      abs (pc.sim$sig - exp.sig) + abs (pc.sim$eps - exp.eps)
    }
    # panchronic
    pc.sim$Pan.dist <- getDist(exp.sig=-1, exp.eps=-1)
    # dead-end
    pc.sim$DE.dist <- getDist(exp.sig=-1, exp.eps=1)
    # ephemeral
    pc.sim$Eph.dist <- getDist(exp.sig=1, exp.eps=1)
    # fuse
    pc.sim$PF.dist <- getDist(exp.sig=1, exp.eps=-1)
    pc.sim
  }
  getGrains <- function (stats, distances, grain=0.1) {
    # create d frame from stats
    d <- data.frame (eps=stats$eps, sig=stats$sig, d=distances)
    # create grains containing coords for each tile
    ps <- seq (-1 + grain, 1, grain)
    grains <- expand.grid (eps=ps, sig=ps)
    # fille tiles with mean value
    for (i in 1:nrow (grains)) {
      higher.eps <- grains[i, 'eps']
      lower.eps <- grains[i, 'eps'] - grain
      eps.pull <- d$eps > lower.eps & d$eps < higher.eps
      higher.sig <- grains[i, 'sig']
      lower.sig <- grains[i, 'sig'] - grain
      sig.pull <- d$sig > lower.sig & d$sig < higher.sig
      grains[i, 'mean'] <- mean (d$d[eps.pull & sig.pull])
      grains[i, 'se'] <- sd (d$d[eps.pull & sig.pull], na.rm = TRUE) /
        sqrt (length (d$d[eps.pull & sig.pull]))
    }
    return (grains)
  }
  pdf (file.path (res.dir, filename), width=9, height=7)
  # remove any that aren't ultrametric or rate.smooted
  if (ignore.chronos) {
    real.stats <- real.stats[real.stats$ultra, ]
  } else {
    real.stats <- real.stats[real.stats$ultra | real.stats$chronos, ]
  }
  # do PCA
  real.stats <- real.stats[!is.na (real.stats$gamma), ]
  real.stats$sig <- NA
  real.stats$eps <- NA
  cols <- c ('sig', 'eps', stat.names)
  input <- rbind (stats[ ,cols], real.stats[, cols])
  pca.res <- prcomp (input[,cols[-c(1,2)]],
                     scale. = TRUE, center = TRUE)
  pca.x.sim <- as.data.frame(pca.res$x[!is.na (input$eps), ])
  pca.x.real <- as.data.frame(pca.res$x[is.na (input$eps), ])
  pca.rot <- as.data.frame (pca.res$rotation)
  prop.var <- round (sapply (pca.res$sdev^2,
                             function (x) Reduce('+', x)/sum (pca.res$sdev^2)), 3)
  names (prop.var) <- colnames (pca.rot)
  # get grains
  pc1.sim <- getGrains (stats, pca.x.sim$PC1, grain=0.1)
  pc2.sim <- getGrains (stats, pca.x.sim$PC2, grain=0.1)
  pc.sim <- data.frame (eps=pc1.sim$eps, sig=pc1.sim$sig,
                        pc1.mean=pc1.sim$mean, pc1.se=pc1.sim$se,
                        pc2.mean=pc2.sim$mean, pc2.se=pc2.sim$se)
  # add real stats
  pc.real <- data.frame (pc1.mean = mean (pca.x.real$PC1),
                         pc2.mean = mean (pca.x.real$PC2),
                         pc1.se = sd (pca.x.real$PC1) /
                           sqrt (length (pca.x.real$PC1)),
                         pc2.se = sd (pca.x.real$PC2) /
                           sqrt (length (pca.x.real$PC2)),
                         Pan.dist=NA, Eph.dist=NA, PF.dist=NA,
                         DE.dist=NA, eps=NA, sig=NA)
  pc.sim <- addDistances(pc.sim)
  res <- rbind (pc.sim, pc.real)
  limitsx <- aes (xmax = pc1.mean + pc1.se,
                  xmin = pc1.mean - pc1.se)
  limitsy <- aes (ymax = pc2.mean + pc2.se,
                  ymin = pc2.mean - pc2.se)
  for (e in c ('Pan.dist', 'Eph.dist', 'PF.dist', 'DE.dist')) {
    p <- ggplot (res, aes_string (x='pc1.mean', y='pc2.mean', colour=e))
    p <- p + geom_point () +
      geom_errorbar(limitsy, width=0.2) +
      geom_errorbarh(limitsx, width=0.2) +
      scale_colour_gradient2 (mid='red', high='blue', na.value='black') +
      xlab (paste0 ('Imbalance - PC1', " (", prop.var['PC1']*100, "%)")) +
      ylab (paste0 ('Loading - PC2', " (", prop.var['PC2']*100, "%)")) +
      theme_bw ()
    print (p)
  }
  closeDevices ()
}

pca <- function (stats, real.stats, stat.names, filename,
                 ignore.chronos=TRUE) {
  pdf (file.path (res.dir, filename), width=9, height=7)
  # remove any that aren't ultrametric or rate.smooted
  if (ignore.chronos) {
    real.stats <- real.stats[real.stats$ultra, ]
  } else {
    real.stats <- real.stats[real.stats$ultra | real.stats$chronos, ]
  }
  real.stats <- real.stats[!is.na (real.stats$gamma), ]
  real.stats$sig <- NA
  real.stats$eps <- NA
  cols <- c ('sig', 'eps', stat.names)
  input <- rbind (stats[ ,cols], real.stats[, cols])
  pca.res <- prcomp (input[,cols[-c(1,2)]],
                     scale. = TRUE, center = TRUE)
  pca.x <- as.data.frame(pca.res$x[!is.na (input$eps), ])
  pca.x$Scenario <- stats$scenario
  pca.x.real <- as.data.frame(pca.res$x[is.na (input$eps), ])
  pca.rot <- as.data.frame (pca.res$rotation)
  prop.var <- round (sapply (pca.res$sdev^2,
                             function (x) Reduce('+', x)/sum (pca.res$sdev^2)), 3)
  names (prop.var) <- colnames (pca.rot)
  # plot means
  res <- data.frame (pca.res$x)
  res$Scenario <- c (stats$scenario, rep ('Emprical', nrow (real.stats)))
  res <- ddply (res, .variables=.(Scenario),
                .fun=summarize, PC1.mean = mean (PC1, na.rm = TRUE),
                PC1.se = sd (PC1, na.rm = TRUE) / sqrt (length (PC1)),
                PC2.mean = mean (PC2, na.rm = TRUE),
                PC2.se = sd (PC2, na.rm = TRUE) / sqrt (length (PC2)))
  limitsx <- aes (colour = Scenario, xmax = PC1.mean + PC1.se,
                 xmin = PC1.mean - PC1.se)
  limitsy <- aes (colour = Scenario, ymax = PC2.mean + PC2.se,
                  ymin = PC2.mean - PC2.se)
  p <- ggplot (res, aes (x=PC1.mean, y=PC2.mean))
  p <- p + geom_point (aes (colour=Scenario)) +
    geom_errorbar(limitsy, width=0.2) +
    geom_errorbarh(limitsx, width=0.2) +
    xlab (paste0 ('Imbalance - PC1', " (", prop.var['PC1']*100, "%)")) +
    ylab (paste0 ('Loading - PC2', " (", prop.var['PC2']*100, "%)")) +
    theme_bw ()
  print (p)
  # plot all points
  comparisons <- list (c ("PC1", "PC2"), c ("PC2", "PC3"), c ("PC1", "PC3"))
  for (comp in comparisons) {
    p <- ggplot (pca.x, aes_string (x = comp[1], y = comp[2])) +
      geom_point (aes (colour=Scenario)) +
      geom_point (data = pca.x.real, colour = 'black', shape = 3) +
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
  closeDevices ()
}

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

tilePlot <- function (stats, distances, grain=0.1, legend.title='d') {
  # convert distances to Z-score
  distances <- (distances - mean (distances)) / sd (distances)
  # create d frame from stats
  d <- data.frame (x=stats$eps, y=stats$sig, d=distances)
  # create p.data containing coords for each tile
  ps <- seq (-1 + grain, 1, grain)
  p.data <- expand.grid (x=ps, y=ps)
  p.data$d <- NA
  # fille tiles with mean value
  for (i in 1:nrow (p.data)) {
    higher.x <- p.data[i, 'x']
    lower.x <- p.data[i, 'x'] - grain
    x.pull <- d$x > lower.x & d$x < higher.x
    higher.y <- p.data[i, 'y']
    lower.y <- p.data[i, 'y'] - grain
    y.pull <- d$y > lower.y & d$y < higher.y
    p.data[i, 'd'] <- mean (d$d[x.pull & y.pull])
  }
  p <- ggplot (p.data, aes (x, y)) + geom_tile (aes (fill=d)) +
    scale_fill_gradient2(low='red', mid='white', high='blue', name=legend.title) +
    labs (x=expression (epsilon), y=expression (sigma)) +
    theme_bw() + theme (axis.title=element_text(size=25))
  return (p)
}

getScenarios <- function (stats) {
  stats$scenario <- NA
  for (i in 1:nrow (stats)) {
    if (stats$eps[i] > 0 & stats$sig[i] > 0) {
      stats$scenario[i] <- 'Eph'
    } else if (stats$eps[i] < 0 & stats$sig[i] > 0) {
      stats$scenario[i] <- 'PF'
    } else if (stats$eps[i] > 0 & stats$sig[i] < 0) {
      stats$scenario[i] <- 'DE'
    } else {
      stats$scenario[i] <- 'Pan'
    }
  }
  return(stats)
}

filter <- function (stats, grain=0.1) {
  max <- 1 - grain
  drop <- NULL
  for (i in 1:nrow (stats)) {
    if (abs (stats$eps[i]) < max | abs (stats$sig[i]) < max) {
      drop <- c (drop, i)
    }
  }
  return(stats[-drop, ])
}

sexdectants <- function (stats) {
  check <- function (min.eps, max.eps, min.sig, max.sig) {
    eps.bool <- stats$eps[i] >= min.eps & stats$eps[i] <= max.eps
    sig.bool <- stats$sig[i] >= min.sig & stats$sig[i] <= max.sig
    return (eps.bool & sig.bool)
  }
  stats$sexdectant <- NA
  for (i in 1:nrow (stats)) {
    if (check (min.eps=0, max.eps=0.5, min.sig=0, max.sig=0.5)) {
      stats$sexdectant[i] <- 'Eph00'
    } else if (check (min.eps=0, max.eps=0.5, min.sig=0.5, max.sig=1)) {
      stats$sexdectant[i] <- 'Eph01'
    } else if (check (min.eps=0.5, max.eps=1, min.sig=0, max.sig=0.5)) {
      stats$sexdectant[i] <- 'Eph10'
    } else if (check (min.eps=0.5, max.eps=1, min.sig=0.5, max.sig=1)) {
      stats$sexdectant[i] <- 'Eph11'
    } else if (check (min.eps=0, max.eps=0.5, min.sig=-1, max.sig=-0.5)) {
      stats$sexdectant[i] <- 'DE00'
    } else if (check (min.eps=0, max.eps=0.5, min.sig=-0.5, max.sig=0)) {
      stats$sexdectant[i] <- 'DE01'
    } else if (check (min.eps=0.5, max.eps=1, min.sig=-1, max.sig=-0.5)) {
      stats$sexdectant[i] <- 'DE10'
    } else if (check (min.eps=0.5, max.eps=1, min.sig=-0.5, max.sig=0)) {
      stats$sexdectant[i] <- 'DE11'
    } else if (check (min.eps=-1, max.eps=-0.5, min.sig=0, max.sig=0.5)) {
      stats$sexdectant[i] <- 'PF00'
    } else if (check (min.eps=-1, max.eps=-0.5, min.sig=0.5, max.sig=1)) {
      stats$sexdectant[i] <- 'PF01'
    } else if (check (min.eps=-0.5, max.eps=0, min.sig=0, max.sig=0.5)) {
      stats$sexdectant[i] <- 'PF10'
    } else if (check (min.eps=-0.5, max.eps=0, min.sig=0.5, max.sig=1)) {
      stats$sexdectant[i] <- 'PF11'
    } else if (check (min.eps=-1, max.eps=-0.5, min.sig=-1, max.sig=-0.5)) {
      stats$sexdectant[i] <- 'Pan00'
    } else if (check (min.eps=-1, max.eps=-0.5, min.sig=-0.5, max.sig=0)) {
      stats$sexdectant[i] <- 'Pan01'
    } else if (check (min.eps=-0.5, max.eps=0, min.sig=-1, max.sig=-0.5)) {
      stats$sexdectant[i] <- 'Pan10'
    } else if (check (min.eps=-0.5, max.eps=0, min.sig=-0.5, max.sig=0)) {
      stats$sexdectant[i] <- 'Pan11'
    } else {
      warning ('Eps and Sig values out of expected range, NAs produced.')
    }
  }
  return(stats)
}

getEDs <- function (trees, scenarios) {
  ed.values <- groups <- NULL
  for (i in 1:length (trees)) {
    res <- calcED(trees[[i]])[ ,1]
    res <- (res - mean (res)) / sd (res)
    ed.values <- append (ed.values, res)
    groups <- append (groups, rep (scenarios[i], length (res)))
  }
  return (data.frame (ed.values, groups))
}

calcTreeStats <- function (trees) {
  engine <- function (i) {
    tree <- trees[[i]]
    # imbalance stats
    ts.tree <- as.treeshape (tree)
    colless.stat <- colless (ts.tree, 'yule')
    sackin.stat <- sackin (ts.tree, 'yule')
    # branching stats
    if (is.null (tree$edge.length)) {
      gamma.stat <- psv.stat <- age <- pd <- NA
    } else {
      # get gamma
      gamma.stat <- gammaStat (tree)
      # to get PSV, create community matrix
      samp <- matrix (rep (1, getSize (trees[[i]]) * 2), nrow=2)
      colnames (samp) <- trees[[i]]$tip.label
      psv.res <- psv (samp, trees[[i]])
      psv.stat <- psv.res[1,1]
      # get age
      age <- getSize (tree, 'rtt')
      # get pd
      pd <- getSize (tree, 'pd')
    }
    data.frame (colless = colless.stat, sackin = sackin.stat,
                gamma = gamma.stat, psv = psv.stat, age, pd)
  }
  mdply (.data = data.frame (i = 1:length (trees)), .fun = engine)[, -1]
}

drawCorresPoints <- function (model, distribution) {
  ## Plot mean, mean and 5% and 95% quantiles of an x distribution's
  ##  equivalent y values
  .draw <- function (actual, predicted, ...) {
    lines (x = c (actual, actual, min.x),
           y = c (min.y, predicted, predicted), ...)
  }
  # find the minimum values of x and y in plotted space
  min.y <- floor (min (model$model[ ,'y'])) + floor (min (model$model[ ,'y']))
  min.x <- floor (min (model$model[ ,'x'])) + floor (min (model$model[ ,'x']))
  # get actual
  q1.x <- quantile (distribution, 0.05, na.rm = TRUE)
  q2.x <- quantile (distribution, 0.95, na.rm = TRUE)
  m.x <- mean (distribution, na.rm = TRUE)
  # get predicted values
  q1.y <- predict (model, data.frame (x = q1.x))
  q2.y <- predict (model, data.frame (x = q2.x))
  m.y <- predict (model, data.frame (x = m.x))
  # draw lines
  .draw (q1.x, q1.y, col = 'red', lty = 2)
  .draw (q2.x, q2.y, col = 'red', lty = 2)
  .draw (m.x, m.y, lwd = 2)
}

calcDistDiff <- function (dist.1, dist.2, n = 1000) {
  # Randomly pull from dists, calc difference, calc mean
  .calc <- function (i) {
    sample.1 <- sample (dist.1, 1)
    sample.2 <- sample (dist.2, 1)
    sample.1 - sample.2
  }
  res <- mdply (.data = data.frame (i = 1:n),
                .fun = .calc)[ ,2]
  mean (res)
}

windowAnalysis <- function (sim.stats, real.stats, size,
                            psi, incr = 0.01,
                            scale = TRUE) {
  .calc <- function (i) {
    samp <- sim.stats[psi < maxs[i] &
                        psi > mins[i]]
    data.frame (diff = abs (calcDistDiff (samp, real.stats)),
                mid = (maxs[i]+mins[i])/2)
  }
  maxs <- seq (from = (min (psi) + size), to = max (psi), by = incr)
  mins <- seq (from = min (psi), to = (max (psi) - size), by = incr)
  res <- mdply (.data = data.frame (i = 1:length (maxs)),
                .fun = .calc)[ ,-1]
  if (scale) {
    res[['diff']] <- res[['diff']] / max (res[['diff']])
  }
  res
}