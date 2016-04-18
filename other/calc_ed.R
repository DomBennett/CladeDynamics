# 04/11/2015
# Calculate rED values

# LIB
library (MoreTreeTools)

# FUNCTIONS
getrEDs <- function (treefiles, treefolder, ...) {
  .get <- function (treefile) {
    .calc <- function (i) {
      EDs <- calcED (trees[[i]])
      rED <- EDs[ ,1] / sum (EDs[ ,1])
      tip <- rownames (EDs)
      data.frame (tip, rED, stringsAsFactors=FALSE)
    }
    tree <- read.tree (file.path (treefolder, treefile))
    if (class (tree) == 'multiPhylo') {
      tree <- tree[[sample (1:length (tree), 1)]]
    }
    if (getSize (tree) > 50) {
      trees <- getSubtrees (tree, min.n=50, max.n=100)
    }
    if (class (trees) == 'multiPhylo') {
      l.data.calc <- data.frame (i=1:length (trees))
      return (mdply (.data=l.data.calc, .fun=.calc)[ ,-1])
    }
  }
  l.data.get <- data.frame (treefile=treefiles)
  mdply (.data=l.data.get, .fun=.get, ...)
}

# EMPIRICAL
treefolder <- file.path ('data', 'parsed_trees_pathD8')
meta <- read.csv (file.path (treefolder, 'treeinfo.csv'), stringsAsFactors=FALSE)
# only pathD8 trees with fewer than 500 tips
treefiles <- meta$filename[meta$rate.smooth == 'pathD8' & meta$ntaxa < 500 & !is.na (meta$ntaxa)]
pathD8.rEDs <- getrEDs (treefiles, treefolder, .progress='time')
# only ultra trees at source
treefiles <- meta$filename[meta$ultra]
ultra.rEDs <- getrEDs (treefiles, treefolder, .progress='time')

# Pan
treefolder <- file.path ('results', 'pan')
treefiles <- list.files (treefolder, pattern='.tre')
pan.rEDs <- getrEDs (treefiles, treefolder, .progress='time')

# DE
treefolder <- file.path ('results', 'de')
treefiles <- list.files (treefolder, pattern='.tre')
de.rEDs <- getrEDs (treefiles, treefolder, .progress='time')

# Hyd
treefolder <- file.path ('results', 'hydra')
treefiles <- list.files (treefolder, pattern='.tre')
hyd.rEDs <- getrEDs (treefiles, treefolder, .progress='time')

# Plot
pathD8.rEDs$Scenario <- 'Emp. (D8)'
ultra.rEDs$Scenario <- 'Emp. (Ultra)'
pan.rEDs$Scenario <- 'Pan'
hyd.rEDs$Scenario <- 'Hyd'
de.rEDs$Scenario <- 'DE'
rEDs <- rbind (pathD8.rEDs, ultra.rEDs, pan.rEDs, hyd.rEDs, de.rEDs)
rEDs$Scenario <- factor (rEDs$Scenario, levels=c ('Emp. (D8)', 'Emp. (Ultra)', 'Pan', 'Hyd', 'DE'))
rEDs$log.rED <- log (rEDs$rED)
p <- ggplot (rEDs, aes (Scenario, log.rED))
p <- p + geom_violin (aes (fill=Scenario)) +
                        theme_bw () +
  theme (axis.text.x = element_blank(), axis.title.x = element_blank(),
         text=element_text(size=25), legend.title=element_blank()) +
  ylab ('log (rED)')
pdf (file.path ('results', 'figures', 'rED.pdf'))
print (p)
dev.off()

# Print max for each
cat ('Max D8 rED: [', max (pathD8.rEDs$rED, na.rm=TRUE), ']\n', sep='')
cat ('Max ultra rED: [', max (ultra.rEDs$rED, na.rm=TRUE), ']\n', sep='')
cat ('Max Pan rED: [', max (pan.rEDs$rED, na.rm=TRUE), ']\n', sep='')
cat ('Max Hyd rED: [', max (hyd.rEDs$rED, na.rm=TRUE), ']\n', sep='')
cat ('Max DE rED: [', max (de.rEDs$rED, na.rm=TRUE), ']\n', sep='')

# Save
save (rEDs, file=file.path ('results', 'rEDs.Rd'))

# load
load (file.path ('results', 'rEDs.Rd'))

names(rEDs)

# Tree size matters
# rEDs$ID <- paste0(rEDs$Scenario, rEDs$treefile, sep='-')
# ntips <- table (rEDs$ID)
# mean.log.rED <- tapply (rEDs$log.rED, rEDs$ID, mean)
# rcols <- rainbow (length (unique (rEDs$Scenario)))
# scenario.counts <- rowSums (table (rEDs$Scenario, rEDs$treefile) > 0)
# cols <- NULL
# for (i in 1:length (rcols)) {
#   cols <- c (cols, rep (rcols[i], scenario.counts[i]))
# }
# plot (as.numeric (mean.log.rED) ~ log(as.numeric(ntips)), col=cols)
# legend ('topright', legend=unique (rEDs$Scenario), col=rcols, pch=19)