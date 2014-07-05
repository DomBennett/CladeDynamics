## 17/06/2014
## D.J. Bennett
## The rise and fall of clades
## Counting success and plotting success of clades

## Deps
library (ape)
library (MoreTreeTools)

## Functions
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

findRiseAndFall <- function (element, min.size = 2) {
  ## Find all clade radiations that rise from and to min.size
  bool <- element >= min.size
  # add shifted bools to get a 1,2 or 3 representing the
  #  position of the number in the sequence
  n.pos <- bool + c (bool[-1], bool[1]) +
    c (bool[length (bool)], bool[-length (bool)])
  # find pos 1s that mark the beginning and end of the radiation
  # with regexp
  n.pos.str <- paste (n.pos, collapse = '')
  start <- regexpr ('123', n.pos.str)
  end <- regexpr ('321', n.pos.str)
  if (any (c (start, end) < 0)) {
    return (NA)
  }
  # return radiation index
  start:(end+2)
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

plotNormalisedSuccess <- function (res, min.time.span = 5, min.size = 2) {
  #   if (is.null (dim (res))) {
  #     stop ('Too few suitable clades -- trying re-running')
  #   }
  n.reps <- ceiling (ncol (res)/3)
  cols <- rep (rainbow (3, alpha = 0.7), n.reps)
  plot (x = 0, y = 0, type = 'n', xlim = c (0.1,1),
        ylim = c (0, 1), main =
          'Normalised Clade Success Through Time',
        xlab = 'Normalised Time', ylab = 'Normalised Success')
  all.x <- all.y <- c ()
  for (i in 1:ncol (res)) {
    # only clades that go from 0 to 0
    # And clades that have minimum time span
    radiation.pos <- findRiseAndFall (res[ ,i], min.size)
    if (is.na (radiation.pos[1]) | length (radiation.pos)
        < min.time.span) {
      next
    }
    # normlise time steps and success
    normalised.time.steps <- 1/length (radiation.pos)
    x <- seq (normalised.time.steps,
              length (radiation.pos)*normalised.time.steps,
              normalised.time.steps)
    y <- res[radiation.pos,i]/max (res[radiation.pos,i])
    lines (x = x, y = y, pch = 19, col = cols[i])
    all.x <- c (all.x, x)
    all.y <- c (all.y, y)
  }
  # add a mean line...
  all.x <- round (all.x, digits = 1)
  uni.all.x <- sort (unique (all.x))
  meaned.y <- rep (NA, length (uni.all.x))
  for (i in 1:length (uni.all.x)) {
    meaned.y[i] <- mean (all.y[all.x %in% uni.all.x[i]])
  }
  lines (x = uni.all.x, y = meaned.y, lwd = 2)
}

plotWithNodeLabels <- function (tree) {
  plot (tree, no.margin = TRUE)
  nodelabels (text = tree$node.label)
}

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
  dev.off()
  system (paste0 ("convert -delay ", delay," *.png ", file.dir))
  file.remove (list.files (pattern = ".png"))
}