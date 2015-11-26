library (MoreTreeTools)
library (grid)
library (gridExtra)

data (EDBMMtrees) # load object
# plot each for comparison
p1<-chromatophylo (edbmm.trees[['DE']][['tree']],
               edge.cols=edbmm.trees[['DE']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
p2<-chromatophylo (edbmm.trees[['Pan']][['tree']],
               edge.cols=edbmm.trees[['Pan']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
p3<-chromatophylo (edbmm.trees[['Hyd']][['tree']],
               edge.cols=edbmm.trees[['Hyd']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
# title settings
tt_theme <- theme(plot.title = element_text(face="bold", hjust=0, vjust=1),
                  text = element_text (size=25))
jpeg (filename="/Users/djb208/Desktop/PRLF_manuscript_paleo/Figure_3.jpg",
     width=1440, height=2880, units="px", pointsize=12,
     quality=100)
grid.arrange(p1 + ggtitle ('A') + tt_theme,
             p2 + ggtitle ('B') + tt_theme,
             p3+ ggtitle ('C') + tt_theme, ncol=1)
dev.off ()

