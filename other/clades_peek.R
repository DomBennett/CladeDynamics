# Check out some of the clades visually
# 30/09/2015

# LIBS
#library (MoreTreeTools)

# PARAMETERS
runlog <- read.csv ('results/clade_results_5/runlog.csv',
                    stringsAsFactors=FALSE)
runlog <- runlog[1:900,]
which (runlog$eps < -0.85 & runlog$sig > 0.85)
treeid <- 401
treeids <- 1:26
pdf ('clade_results.pdf')
for (treeid in treeids) {
  # INPUT
  c.file <- paste0 ("results/clade_results_5/tree", treeid, "_clades_ind.csv")
  cs.file <- paste0 ("results/clade_results_5/tree", treeid, "_clade_ind_stats.csv")
  trees.file <- paste0 ("results/clade_results_5/tree", treeid, "_recorded_trees.rda")
  clades <- read.csv (c.file)
  clade.stats <- read.csv (cs.file)
  #load (trees)
  
  # LOOK-SEE
  selection <- clade.stats$cid [clade.stats$end < max (clade.stats$end) &
                                  clade.stats$max.size > 5 &
                                  clade.stats$start > 25]
  #plotClades (clades, cids=selection)
  p <- plotClades (clades, cids=selection, merge=TRUE)
  sig <- signif (runlog[treeid,'sig'], 2)
  eps <- signif (runlog[treeid,'eps'], 2)
  p <- p + ggtitle (paste0 ('Sig = ', sig, ' | Eps = ', eps))
  print (p)
}
dev.off()

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
