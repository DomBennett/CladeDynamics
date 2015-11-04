# 04/11/2015
# Calculate rED values

# LIB
library (MoreTreeTools)

# FUNCTIONS
getrEDs <- function (treefiles, treefolder, ...) {
  .calc <- function (treefile) {
    tree <- read.tree (file.path (treefolder, treefile))
    if (class (tree) == 'multiPhylo') {
      tree <- tree[[sample (1:length (tree), 1)]]
    }
    EDs <- calcED (tree)
    rED <- EDs[ ,1] / sum (EDs[ ,1])
    tip <- rownames (EDs)
    data.frame (tip, rED, stringsAsFactors=FALSE)
  }
  l.data <- data.frame (treefile=treefiles)
  mdply (.data=l.data, .fun=.calc, ...)
}

# EMPIRICAL
treefolder <- file.path ('data', 'parsed_trees_pathD8')
meta <- read.csv (file.path (treefolder, 'treeinfo.csv'), stringsAsFactors=FALSE)
# only pathD8 trees with fewer than 500 tips
treefiles <- meta$filename[meta$rate.smooth == 'pathD8' & meta$ntaxa < 500 & !is.na (meta$ntaxa)]
emp.rEDs <- getrEDs (treefiles, treefolder, .progress='time')

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
emp.rEDs$Scenario <- 'Emp.'
pan.rEDs$Scenario <- 'Pan'
hyd.rEDs$Scenario <- 'Hyd'
de.rEDs$Scenario <- 'DE'
rEDs <- rbind (emp.rEDs, pan.rEDs, hyd.rEDs, de.rEDs)
rEDs$Scenario <- factor (rEDs$Scenario, levels=c ('Emp.', 'Pan', 'Hyd', 'DE'))
p <- ggplot (rEDs, aes (Scenario, rED))
p <- p + geom_violin (aes (fill=Scenario)) +
                        theme_bw () +
  theme (axis.text.x = element_blank(), axis.title.x = element_blank(),
         text=element_text(size=25), legend.title=element_blank())
pdf (file.path ('results', 'figures', 'rED.pdf'))
print (p)
dev.off()

# Print max for each
cat ('Max Emp rED: [', max (emp.rEDs$rED, na.rm=TRUE), ']\n', sep='')
cat ('Max Pan rED: [', max (pan.rEDs$rED, na.rm=TRUE), ']\n', sep='')
cat ('Max Hyd rED: [', max (hyd.rEDs$rED, na.rm=TRUE), ']\n', sep='')
cat ('Max DE rED: [', max (de.rEDs$rED, na.rm=TRUE), ']\n', sep='')

# Save
save (rEDs, file=file.path ('results', 'rEDs.Rd'))
