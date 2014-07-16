## 17/06/2014
## D.J. Bennett
## The rise and fall of clades
## Counting success and plotting success of clades

## Deps
library (ape)
library (MoreTreeTools)
library (ggplot2)

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

plotNormalisedSuccess <- function (res, min.time.span = 5, min.size = 2) {
  #   if (is.null (dim (res))) {
  #     stop ('Too few suitable clades -- trying re-running')
  #   }
  n.reps <- ceiling (ncol (res)/3)
  cols <- rep (rainbow (3, alpha = 0.7), n.reps)
  plot (x = 0, y = 0, type = 'n', xlim = c (-1,1),
        ylim = c (0, 1), main =
          'Normalised Clade Success Through Time',
        xlab = 'Normalised Time', ylab = 'Normalised Success')
  all.x <- all.y <- c ()
  for (i in 1:ncol (res)) {
    normalised.element <- normalise (res[ ,i], min.size, min.time.span)
    if (any (is.na (normalised.element))) {
      next
    }
    lines (x = normalised.element$x, y = normalised.element$y,
           pch = 19, col = cols[i])
    all.x <- c (all.x, normalised.element$x)
    all.y <- c (all.y, normalised.element$y)
  }
  if (any (is.null (all.x))) {
    print ('No clades met criteria minimum criteria for normalised plotting.')
  } else {
    # add a mean line...
    all.x <- round (all.x, digits = 1)
    uni.all.x <- sort (unique (all.x))
    meaned.y <- rep (NA, length (uni.all.x))
    for (i in 1:length (uni.all.x)) {
      meaned.y[i] <- mean (all.y[all.x %in% uni.all.x[i]])
    }
    lines (x = uni.all.x, y = meaned.y, lwd = 2)
  }
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
    res <- calcFairProportion (trees[[i]])
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

histClageAges <- function (res.dirs) {
  ## Take a vector of res dirs and plot clade ages in hist
  ages <- bias <- c ()
  for (each in res.dirs) {
    res <- read.csv (file = file.path (
      'results', each, 'clades_through_time.csv'))
    ages <- c (ages, colSums (res != 0))
    each.bias <- substr (x = each,
                         start = (regexpr ("_bias", each) + 5),
                         stop = (regexpr ("_date", each) - 1))
    bias <- c (bias, rep (each.bias, ncol (res)))
  }
  clade.ages <- data.frame (ages = ages, bias = bias)
  hist <- ggplot(clade.ages, aes (x = ages, fill = bias))
  hist + geom_bar()
}