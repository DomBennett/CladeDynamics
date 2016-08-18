library (MoreTreeTools)
library (grid)
library (gridExtra)


# Function from MoreTreeTools
chromatophylo2 <- function (tree, edge.cols = NULL, legend.title = "", 
          reduce.overlap = TRUE) 
{
  if (!is.binary.tree(tree)) {
    stop("Tree must be bifurcating")
  }
  N <- length(tree$tip.label)
  tree.age <- getSize(tree, "rtt")
  x.spacer <- tree.age/100
  y.spacer <- 1/N
  x1 <- tree.age
  node <- length(tree$tip.label) + 1
  edge <- which(tree$edge[, 1] == node)
  p.env <- new.env(parent = emptyenv())
  p.env$N <- N
  p.env$all.edges <- 1:nrow(tree$edge)
  p.env$tree.stats <- getTreeStats(tree)
  p.env$p.data <- data.frame(x = c(x1, x1), y = c(0, 1), edge = c(0, 
                                                                  0))
  p.env$y.spacer <- y.spacer
  p.env$x.spacer <- x.spacer
  p.env$tree <- tree
  p.env$root.node <- length(tree$tip.label) + 1
  MoreTreeTools:::.cpMkPData(x1, c(0, 1), edge, p.env)
  p.env$p.data$line <- rep(1:(nrow(p.env$p.data)/2), each = 2)
  p.env$lines <- unique(p.env$p.data$line)[-1]
  edges <- p.env$p.data$edge[match(p.env$lines, p.env$p.data$line)]
  p.env$maxmin <- data.frame(line = p.env$lines, edge = edges, 
                             max.x = NA, min.x = NA, max.y = NA, min.y = NA)
  MoreTreeTools:::.cpGetMaxMin(p.env$lines, p.env)
  while (reduce.overlap) {
    p.env$overlap <- FALSE
    plyr::m_ply(.data = data.frame(l = p.env$lines), .fun = .cpCheckLine, 
                p.env = p.env)
    if (!p.env$overlap) {
      break
    }
  }
  p.data <- p.env$p.data
  if (!is.null(edge.cols)) {
    matched.is <- match(p.data$edge, edge.cols$edge)
    p.data$col <- edge.cols$col[matched.is]
  }
  p <- ggplot(p.data, aes(x = x, y = y, group = edge))
  p <- p + geom_line(aes(colour = col), size = 0.25) + 
    scale_colour_gradient2(low="white", high="red")
  theme.opts <- theme(panel.background = element_rect(fill = "gray15"), 
                      axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
                      axis.text.y = element_blank(), panel.grid.major.y = element_blank(), 
                      panel.grid.minor.y = element_blank(), panel.grid.major.x = element_line(colour = "black"), 
                      panel.grid.minor.x = element_line(colour = "gray5"))
  p <- p + scale_x_reverse() + xlab("Time") + theme.opts
  p
}

data (EDBMMtrees) # load object from MoreTreeTools
# https://github.com/DomBennett/MoreTreeTools/blob/master/other/model_edbmm.R
# plot each for comparison
de_p<-chromatophylo2(tree=edbmm.trees[['DE']][['tree']],
               edge.cols=edbmm.trees[['DE']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
pn_p<-chromatophylo2 (edbmm.trees[['Pan']][['tree']],
               edge.cols=edbmm.trees[['Pan']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
hyd_p<-chromatophylo2 (edbmm.trees[['Hyd']][['tree']],
               edge.cols=edbmm.trees[['Hyd']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
# title settings
tt_theme <- theme(plot.title=element_text(hjust=0, vjust=1),
  legend.position='none', axis.title.x=element_blank(),
  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
  text = element_text(size=8))
tiff("~/Desktop/figure.tiff", width=18, height=9, units="cm",
     res=1200)
grid.arrange(de_p + ggtitle ('(DE)') + tt_theme,
             hyd_p + ggtitle ('(Hydra)') + tt_theme,
             pn_p + ggtitle ('(Pan)') + tt_theme, ncol=3)
dev.off ()
