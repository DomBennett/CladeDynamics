# Read in latest results for clade analysis for DE, Pan and Hydra scenarios

library (grid)
library (gridExtra)
library(MoreTreeTools)
library(Hmisc)

plotClades2 <- function (clades, N = 3, clade.stats = NULL, cids = NULL, legend = FALSE, 
          merge = FALSE, lty=1) 
{
  .combine <- function(cid) {
    n <- as.vector(clades[, cid])
    t <- 1:length(n)
    data.frame(n, t, c = cid)
  }
  .normalise <- function(cid) {
    n <- as.vector(clades[, cid])
    n <- n[n != 0]
    n <- n/min(n)
    n <- n/max(n)
    t <- seq(0, 1, length.out = length(n))
    data.frame(n, t)
  }
  if (is.null(cids)) {
    if (merge) {
      cids <- colnames(clades)
    }
    else {
      if (is.null(clade.stats)) {
        clade.stats <- calcCladeStats(clades)
      }
      cids <- order(clade.stats$tot.size, decreasing = TRUE)[1:N]
      cids <- which(clade.stats$cid %in% colnames(clades))
    }
  }
  gt <- paste0("N = ", length(cids))
  p.data <- data.frame(cid = cids, stringsAsFactors = FALSE)
  if (merge) {
    p.data <- plyr::mdply(.data = data.frame(p.data), .fun = .normalise)[, 
                                                                         -1]
    dps <- round(log(x = nrow(clades), base = 10))
    p.data$t <- signif(p.data$t, dps)
    p.data <- plyr::ddply(p.data, plyr::.(t), plyr::summarise, 
                          mean.n = mean(n), upper.sd = mean(n) + sd(n), lower.sd = mean(n) - 
                            sd(n))
    p <- ggplot(p.data) + geom_line(aes(x = t, y = mean.n), 
                                    colour = "black") +
      geom_line(aes(x = t, y = upper.sd), 
                colour = "red", lty = 3) +
      geom_line(aes(x = t, y = lower.sd), 
                colour = "red", lty = 3) + theme_bw() + xlab("Normalised time step") + 
      ylab("Normalised N")
  }
  else {
    p.data <- plyr::mdply(.data = p.data, .fun = .combine)[, -1]
    p <- ggplot(p.data, aes(x = t, y = n, colour = c)) + 
      geom_line(size=lty) + theme_bw() + xlab("Time step") + ylab("N")
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p + ggtitle(gt))
}


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

# plot examples of each scanerio
ps <- list()
scenarios <- c("Pan", "Hydra", "DE")
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
      p <- plotClades2(clades, cids=cids[1:100], merge=FALSE, lty=.1)
      p <- p + ggtitle(paste0("(", capitalize(s), ")"))
      p <- p + xlab("") + ylab("") + ylim(0, 25)
      ps <- c(list(p), ps)
      break
    }
  }
}

tt <- theme(text=element_text(size=8, hjust=0), plot.margin=unit(c(.1,.2,-.1,-.1), "cm"))
ttt <- theme(axis.text.y=element_blank(), text=element_text(size=8, hjust=0),
             plot.margin=unit(c(.1,.2,-.1,-.3), "cm"))
tiff("~/Desktop/figure.tiff", width=18, height=9, units="cm",
     res=1200)
grid.arrange(ps[[1]] + tt,
             ps[[2]] + ttt,
             ps[[3]] + ttt,
             ncol=3)
dev.off()
