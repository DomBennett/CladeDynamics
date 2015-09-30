# Check out some of the clades visually
# 30/09/2015

# LIBS
library (MoreTreeTools)

# PARAMETERS
treeid <- 5

# INPUT
runlog <- read.csv ('results/analysis_5/runlog.csv',
                    stringsAsFactors=FALSE)
c.file <- paste0 ("results/analysis_5/tree", treeid, "_clades.csv")
cs.file <- paste0 ("results/analysis_5/tree", treeid, "_clade_stats.csv")
trees.file <- paste0 ("results/analysis_5/tree", treeid, "_recorded_trees.rda")
clades <- read.csv (c.file)
clade.stats <- read.csv (cs.file)
#load (trees)

# LOOK-SEE
selection <- clade.stats$cid [clade.stats$end < max (clade.stats$end) &
                                clade.stats$max.size > 5 &
                                clade.stats$start > 25]
#plotClades (clades, cids=selection)
plotClades (clades, cids=selection, merge=TRUE)
print (runlog[treeid,])
