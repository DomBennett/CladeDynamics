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