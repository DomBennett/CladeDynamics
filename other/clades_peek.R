# Check out some of the clades visually
# 30/09/2015

# LIBS
library (MoreTreeTools)

# PARAMETERS
runlog <- read.csv ('results/clade_results_5/runlog.csv',
                    stringsAsFactors=FALSE)
runlog <- runlog[1:900,]
# get pand, de and h scernarios
treeids <- which (runlog$sig < -0.8 & runlog$eps < -0.8)
#treeids <- c(treeids, which (runlog$sig < -0.8 & runlog$eps > 0.8))
#treeids <- c(treeids, which (runlog$sig < -0.8 & runlog$eps < 0.1 & runlog$eps > -0.1))
for (treeid in treeids[2:length(treeids)]) {
  # INPUT
  c.file <- paste0 ("results/clade_results_5/tree", treeid, "_clades_ind.csv")
  cs.file <- paste0 ("results/clade_results_5/tree", treeid, "_clade_ind_stats.csv")
  trees.file <- paste0 ("results/clade_results_5/tree", treeid, "_recorded_trees.rda")
  tmp_clds <- read.csv (c.file, stringsAsFactors=FALSE)
  colnames(tmp_clds) <- paste0(treeid, '_', colnames(tmp_clds))
  tmp_cldstts <- read.csv (cs.file, stringsAsFactors=FALSE)
  tmp_cldstts$cid <- paste0(treeid, '_', tmp_cldstts$cid)
  clades <- cbind(clades, tmp_clds)
  clade.stats <- rbind(clade.stats, tmp_cldstts)
}
pdf ('clade_results.pdf')
cids <- clade.stats$cid[clade.stats$max.size > 15 &
                          clade.stats$start > 25 &
                          clade.stats$end < 75]
if(length(cids) > 1) {
  p <- plotClades (clades, cids=cids, merge=TRUE, nbin=10)
  sig <- signif (runlog[treeid,'sig'], 2)
  eps <- signif (runlog[treeid,'eps'], 2)
  p <- p + ggtitle (paste0 ('Sig = ', sig, ' | Eps = ', eps))
  print (p)
}
dev.off()


plotClades <- function (clades, N=3, clade.stats=NULL, cids=NULL,
                        legend=FALSE, merge=FALSE, nbin=5) {
    # Internals
    .bin <- function(n, t, ints, .fun){
      .tmean <- function(t1, t2) {
        bool <- t < t2 & t >= t1
        .fun(n[bool], na.rm=TRUE)
      }
      l_data <- data.frame(t1=ints[1:(length(ints) - 1)],
                           t2=ints[2:length(ints)])
      mdply(l_data, .fun=.tmean)
    }
    .combine <- function (cid) {
      # combine into single dataframe for ggplot
      n <- as.vector (clades[ ,cid])
      t <- 1:length(n)
      ints <- seq(1, length(n), length.out=nbin+1)
      binned.n <- .bin(n, t, ints, .fun=mean)
      t <- binned.n$t1 + (binned.n$t2 - binned.n$t1)/2
      n <- binned.n[ ,'V1']
      data.frame (n, t, c=cid)
    }
    .normalise <- function (cid) {
      n <- as.vector (clades[ ,cid])
      n <- n[n != 0]
      n <- n / min (n)
      n <- n / max(n)
      t <- seq(0, 1, length.out=length(n))
      data.frame (n, t)
    }
    if (is.null (cids)) {
      if (merge) {
        # if no cids given and merge, use all clades
        cids <- colnames (clades)
      } else {
        # otherwise use top N clades
        if (is.null (clade.stats)) {
          clade.stats <- calcCladeStats (clades)
        }
        # get top N clades by total size
        cids <- order (clade.stats$tot.size, decreasing=TRUE)[1:N]
        # get their position in clades
        cids <- which (clade.stats$name[i] %in% colnames (clades))
      }
    }
    gt <- paste0 ('N = ', length (cids))
    p.data <- data.frame(cid=cids, stringsAsFactors=FALSE)  # avoid levels
    if (merge) {
      # if merged into single plot
      n.data <- mdply(.data=data.frame(p.data), .fun=.normalise)[ ,-1]
      ints <- seq(0, 1, length.out=nbin + 1)
      means <- .bin(n.data$n, t=n.data$t, ints=ints, .fun=mean)
      sds <- .bin(n.data$n, t=n.data$t, ints=ints, .fun=sd)
      t <- means$t1 + (means$t2 - means$t1)/2
      p.data <- data.frame(t=t)
      p.data$mean.n <- means[ ,'V1']
      p.data$upper.sd <- means[ ,'V1'] + sds[ ,'V1']
      p.data$lower.sd <- means[ ,'V1'] - sds[ ,'V1']
      p <- ggplot (p.data) +
        geom_line (aes (x=t, y=mean.n), colour='black') +
        geom_line (aes (x=t, y=upper.sd), colour='red', lty=3) +
        geom_line (aes (x=t, y=lower.sd), colour='red', lty=3) +
        theme_bw() + xlab ('Normalised time step') + ylab ('Normalised N')
    } else {
      # combine clades into ggplot plot dataframe
      p.data <- mdply (.data=p.data, .fun=.combine)[ ,-1]
      p <- ggplot (p.data, aes (x=t, y=n, colour=c)) +
        geom_line() + theme_bw() +
        xlab ('Time step') + ylab ('N')
    }
    if (!legend) {
      p <- p + theme (legend.position="none")
    }
    return (p + ggtitle(gt))
}

p.pan <- p
p.de <- p

pdf ('example_clade_plots.pdf')
print (p.pan + ggtitle ('Pan'))
print (p.de + ggtitle ('DE'))
dev.off()

# using CM and CG
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'clade_tools.R'))
stats <- readIn ('analysis_5')
stats <- stats[1:1000,]
clade.stats <- readInCladeStats (stats=stats, 'clade_results_5')
pull <- !is.na (clade.stats$cm)
p <- tilePlot (clade.stats[pull,], clade.stats$cm[pull], legend.title='CM, Z-score',
               grain=0.25)
print (p)
pull <- !is.na (clade.stats$cg)
p <- tilePlot (clade.stats[pull,], clade.stats$cg[pull], legend.title='CG, Z-score')
print (p)
mean (clade.stats$cg[pull])
# base stats
tapply (clade.stats$cm, clade.stats$scenario, mean)
tapply (clade.stats$cg, clade.stats$scenario, mean)
tapply (clade.stats$max.size, clade.stats$scenario, mean)
tapply (clade.stats$max.size, clade.stats$scenario, sd)
tapply (clade.stats$time.span, clade.stats$scenario, mean)
tapply (clade.stats$time.span, clade.stats$scenario, sd)
