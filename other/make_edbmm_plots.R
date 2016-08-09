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
tt_theme <- theme(plot.title=element_text(hjust=0, vjust=1),
  legend.position='none', axis.title.x=element_blank(),
  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
  text = element_text (size=25))
tiff(filename="/Users/djb208/Desktop/Figure_5.tiff",
     width=1440, height=2880, units="px", pointsize=12)
grid.arrange(p1 + ggtitle ('DE') + tt_theme,
             p2 + ggtitle ('Hydra') + tt_theme,
             p3 + ggtitle ('Pan') + tt_theme, ncol=1)
dev.off ()

