library (MoreTreeTools)
library (grid)
library (gridExtra)

data (EDBMMtrees) # load object
# plot each for comparison
de_p<-chromatophylo (edbmm.trees[['DE']][['tree']],
               edge.cols=edbmm.trees[['DE']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
pn_p<-chromatophylo (edbmm.trees[['Pan']][['tree']],
               edge.cols=edbmm.trees[['Pan']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
hyd_p<-chromatophylo (edbmm.trees[['Hyd']][['tree']],
               edge.cols=edbmm.trees[['Hyd']][['edge.diversity']],
               legend.title='Diversity', reduce.overlap=FALSE)
# title settings
tt_theme <- theme(plot.title=element_text(hjust=0, vjust=1),
  legend.position='none', axis.title.x=element_blank(),
  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
  text = element_text (size=25))
nw_col <- scale_colour_gradient2(low="white", high="red")
tiff(filename="/Users/djb208/Desktop/Figure_5.tiff",
     width=2880, height=1440, units="px", pointsize=12)
grid.arrange(de_p + nw_col + ggtitle ('(DE)') + tt_theme,
             hyd_p + nw_col + ggtitle ('(Hydra)') + tt_theme,
             pn_p + nw_col + ggtitle ('(Pan)') + tt_theme, ncol=3)
dev.off ()
