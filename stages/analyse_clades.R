# Read in latest results for clade analysis for DE, Pan and Hydra scenarios

library(MoreTreeTools)

all_clade_stats <- data.frame(cid=NA, tot.size=NA, max.size=NA,
                              time.span=NA, start=NA, end=NA,
                              cm=NA, cg=NA, scenario=NA)
scenarios <- c('de', 'pan', 'hydra')
for(s in scenarios) {
  cat('... ', s, '\n')
  input_dir <- file.path('results', 'clade_res', s, 'clade_results')
  fls <- list.files(input_dir)
  fls <- fls[grepl('*.rda', fls)]
  treeids <- gsub('_.*', '', fls)
  for(treeid in treeids) {
    clade_stats <- read.csv(file.path(input_dir, paste0(treeid, '_clade_stats.csv')),
                            stringsAsFactors=FALSE)[ ,-1]
    clade_stats$scenario <- s
    all_clade_stats <- rbind(all_clade_stats, clade_stats)
  }
}

# CALC CENTRE OF MASS
all_clade_stats <- all_clade_stats[-1, ]
bool <- all_clade_stats$tot.size > 100 &
  all_clade_stats$start > 25 &
  all_clade_stats$end < 75
sum(bool)
table(all_clade_stats$scenario)
tapply(all_clade_stats$cm,
       all_clade_stats$scenario,
       mean, na.rm=TRUE)
tapply(all_clade_stats$cm[bool],
       all_clade_stats$scenario[bool],
       sd, na.rm=TRUE)
tapply(all_clade_stats$cm[bool],
       all_clade_stats$scenario[bool],
       sd, na.rm=TRUE)
tapply(all_clade_stats$time.span[bool],
       all_clade_stats$scenario[bool],
       mean, na.rm=TRUE)
tapply(all_clade_stats$time.span[bool],
       all_clade_stats$scenario[bool],
       sd, na.rm=TRUE)
tapply(all_clade_stats$max.size[bool],
       all_clade_stats$scenario[bool],
       mean, na.rm=TRUE)
tapply(all_clade_stats$max.size[bool],
       all_clade_stats$scenario[bool],
       sd, na.rm=TRUE)
tapply(all_clade_stats$tot.size[bool],
       all_clade_stats$scenario[bool],
       mean, na.rm=TRUE)
tapply(all_clade_stats$tot.size[bool],
       all_clade_stats$scenario[bool],
       sd, na.rm=TRUE)
tapply(all_clade_stats$cg[bool],
       all_clade_stats$scenario[bool],
       mean, na.rm=TRUE)
tapply(all_clade_stats$cg[bool],
       all_clade_stats$scenario[bool],
       sd, na.rm=TRUE)
save(all_clade_stats, file='all_clade_stats.RData')
# load('all_clade_stats.RData')


pdf('clade_res.pdf')
for(s in scenarios) {
  input_dir <- file.path('results', 'clade_res', s, 'clade_results')
  fls <- list.files(input_dir)
  fls <- fls[grepl('*.rda', fls)]
  treeids <- gsub('_.*', '', fls)
  treeids <- sample(treeids)
  for(treeid in treeids) {
    clade_stats <- read.csv(file.path(input_dir, paste0(treeid, '_clade_stats.csv')),
                            stringsAsFactors=FALSE)
    print(max(clade_stats$tot.size))
    if(max(clade_stats$tot.size) > 10000) {
      clades <- read.csv(file.path(input_dir, paste0(treeid, '_clades.csv')),
                         stringsAsFactors=FALSE)
      ord <- order(clade_stats$tot.size, decreasing=TRUE)
      clade_stats <- clade_stats[ord, ]
      cids <- clade_stats$cid[clade_stats$max.size < 25 &
                                clade_stats$start > 25 &
                                clade_stats$end < 75]
      p <- plotClades(clades, cids=cids[1:100], merge=FALSE, nbin=10)
      p <- p + ggtitle(s)
      print (p)
      break
    }
  }
}
dev.off()