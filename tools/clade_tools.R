# load compare_tools first

readInCladeStats <- function (stats, analysis.name,
                              filter=TRUE, min.size=10) {
  # search results path for Rd, else read in each csv
  rdfile <- file.path ('results', analysis.name, 'clade_stats.Rd')
  if (file.exists (rdfile)) {
    load (rdfile)
    return (clade.stats)
  }
  # plyr loop
  .run <- function (i) {
    if (file.exists (file.path ('results', analysis.name,
                                 cladestats.files[i]))) {
      part.clade.stats <- read.csv (file.path ('results', analysis.name,
                                               cladestats.files[i]))[,-1]
      part.clade.stats$sig <- stats$sig[i]
      part.clade.stats$eps <- stats$eps[i]
      if (filter) {
        # ... ignore clades that started
        part.clade.stats <- part.clade.stats[part.clade.stats$start > 25,]
        # ... ignore clades that were still extant
        part.clade.stats <- part.clade.stats[part.clade.stats$end != max (part.clade.stats$end),]
        # ... only take clades over a certain size
        part.clade.stats <- part.clade.stats[part.clade.stats$tot.size > min.size, ]
      }
      part.clade.stats <- getScenarios (part.clade.stats)
      # bind+push
      clade.stats <<- rbind (clade.stats, part.clade.stats)
    }
  }
  # get all files
  cladestats.files <- sub ('\\.tre', '_clade_ind_stats.csv',
                           stats$treefilename)
  # template
  clade.stats <- data.frame (cid=NA, tot.size=NA,
                             max.size=NA, time.span=NA,
                             start=NA, end=NA, cm=NA, cg=NA,
                             sig=NA, eps=NA, scenario=NA)
  m_ply (.data=data.frame (i=1:nrow (stats)),
         .fun=.run)
  clade.stats <- clade.stats[-1, ]
  save (clade.stats, file=rdfile)
  return (clade.stats)
}

## Clade analysis functions with tests (no longer part of pipeline 05/09/2014)
findNonZeros <- function (element) {
  ## Find all results that are not zero
  ##  save for those that occur 1 before
  ##  and 1 after
  bool <- element > 0
  # add shifted bools to get a 1,2 or 3 representing the
  #  position of the number in the sequence
  n.pos <- bool + c (bool[-1], bool[1]) +
    c (bool[length (bool)], bool[-length (bool)])
  # return a bool of all n.pos above 0
  n.pos > 0
}

plotSuccess <- function (res, section = 'all',
                         thinning = 'all') {
  # Plot as separate lines for the success of each
  # clade through time
  # section = 'all' -- for all clades
  # else do for bottom number
  # thinning allows only a proportion of the
  # clade results to be plotted
  if (section != 'all') {
    top <- ncol (res)
    bottom <- top - section
    res <- res[bottom:top]
  }
  if (thinning != 'all') {
    # choose those to plot by biasing towards the
    # older ones -- makes it easier to see
    # too many young ones
    random.is <- sample (1: ncol (res),
                         thinning, prob = 1/1: ncol (res))
    res <- res[ ,random.is]
  }
  n.clades <- ncol (res)
  # re-use the first 3 rainbow colours
  # red, green and blue
  n.reps <- ceiling (n.clades/3)
  cols <- rep (rainbow (3, alpha = 0.7), n.reps)
  x.max <- max (as.numeric (rownames (res)))
  x.min <- min (as.numeric (rownames (res)))
  y.max <- max (res)
  y.min <- min (res)
  plot (0, type = 'n', xlim = c (x.min, x.max), ylim =
          c (y.min, y.max), ylab = 'N', xlab = 'Time',
        main = 'Clade success through time')
  for (i in 1:n.clades) {
    # don't plot all 0s -- otherwise very cluttered
    # but do plot the zero for start and end of clade
    non.zero <- findNonZeros (res[ ,i])
    y <- res[non.zero, i]
    x <- as.numeric (rownames (res)[non.zero])
    lines (y = y, x = x, col = cols[i], pch = 19)
  }
}

normalise <- function (element, min.size, min.time) {
  ## Return a list of normalised clade success and time
  ##  centred at its peak
  if (element[length (element)] != 0 |
        element[1] != 0) { # only extinct clades
    return (NA)
  }
  if (max (element) < min.size) {
    return (NA)
  }
  if (sum (element != 0) < min.time) {
    return (NA)
  }
  # remove excessive zeros
  y <- element[findNonZeros (element)]
  # find peak(s)
  peak <- which (y == max (y))
  peak <- round (mean (peak))
  # centre time around peak
  x <- seq ((0 - (peak-1)), (length (y) - peak))
  # normalise
  x <- x/max (x)
  y <- y/max (y)
  list (x = x, y = y)
}

plotNormalisedSuccess <- function (res, min.time.span=5, min.size=2,
                                   plot.all=FALSE, plot.sd=FALSE,
                                   new.plot=TRUE, col='black') {
  n.reps <- ceiling (ncol (res)/3)
  cols <- rep (rainbow (3, alpha = 0.7), n.reps)
  if (new.plot) {
    plot (x = 0, y = 0, type = 'n', xlim = c (-1,1),
          ylim = c (0, 1), main =
            'Normalised Clade Success Through Time',
          xlab = 'Normalised Time', ylab = 'Normalised Success')
  }
  all.x <- all.y <- c ()
  for (i in 1:ncol (res)) {
    normalised.element <- normalise (res[ ,i], min.size, min.time.span)
    if (any (is.na (normalised.element))) {
      next
    }
    if (plot.all) {
      lines (x = normalised.element$x, y = normalised.element$y,
             pch = 19, col = cols[i])
    }
    all.x <- c (all.x, normalised.element$x)
    all.y <- c (all.y, normalised.element$y)
  }
  if (any (is.null (all.x))) {
    cat ('\n... No clades met minimum criteria for normalised plotting')
  } else {
    # add a mean line...
    all.x <- round (all.x, digits = 1)
    uni.all.x <- sort (unique (all.x))
    meaned.y <- sd.y <- rep (NA, length (uni.all.x))
    for (i in 1:length (uni.all.x)) {
      meaned.y[i] <- mean (all.y[all.x %in% uni.all.x[i]])
      sd.y[i] <- sd (all.y[all.x %in% uni.all.x[i]])
    }
    lines (x = uni.all.x, y = meaned.y, lwd = 2, col=col)
    if (plot.sd) {
      lines (x = uni.all.x, y = meaned.y + sd.y,
             col='red', lty=3)
      lines (x = uni.all.x, y = meaned.y - sd.y,
             col='red', lty=3)
    }
  }
}