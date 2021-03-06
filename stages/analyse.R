# 16/07/2014
# D.J. Bennett
# Comparing simulated and empirical trees -- Interactive script

# LIBS
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))
library (gridExtra)
library (outliers)

# DIRS
data.dir <- file.path ('data', 'treestats')

# PARAMETERS
analysis.name <- 'analysis_5'
res.dir <- file.path ('results', 'figures')
data.dir <- file.path ('data', 'treestats')
empirical.filenames <- c ('min50_max500_rspathD8.Rd',
                          'min50_max500_rschronopl.Rd',
                          'min50_max500_rschronoMPL.Rd')

# INPUT
stats <- readIn (analysis.name)
extreme <- rbind (readIn ('Pan'), readIn ('Eph'),
                  readIn ('DE'), readIn ('PF'), readIn ('Hydra'))
# load pre-calculated empirical tree stats -- real.stats
real.stats <- readInMultiple (file.path (data.dir, empirical.filenames))
# add taxoinfo to real.stats
taxoinfo <- read.csv (file.path (data.dir, 'taxoinfo.csv'),
                      stringsAsFactors=FALSE)
real.stats$phylum <- real.stats$class <- real.stats$order <- NA
taxoinfo <- taxoinfo[taxoinfo$treefile %in% real.stats$filename, ]
match.is <- match (taxoinfo$treefile, real.stats$filename)
real.stats$phylum[match.is] <- taxoinfo$phylum
real.stats$class[match.is] <- taxoinfo$class
real.stats$order[match.is] <- taxoinfo$order
real.stats$phylum[is.na (real.stats$phylum)] <- 'Unknown'
real.stats$class[is.na (real.stats$class)] <- 'Unknown'
real.stats$order[is.na (real.stats$order)] <- 'Unknown'

names(extreme)

# PARSE
cat ('\nDropping outliers ....')
# remove branch results where psv is greater than 1 -- this is impossible!
# also drop the gamma for the same tree
psvs <- which (grepl ('psv', colnames (real.stats)))
gammas <- which (grepl ('gamma', colnames (real.stats)))
for (i in 1:length (psvs)) {
  psv.i <- psvs[i]
  gamma.i <- gammas[i]
  pull <- real.stats[ ,psv.i] > 1 & !is.na (real.stats[ ,psv.i])
  cat (paste0 (colnames (real.stats)[psv.i], " -- ", sum (pull), '\n'))
  real.stats[pull, psv.i] <- NA
  real.stats[pull, gamma.i] <- NA
}
# Extremely conservative removal of outliers
real.stats <- dropOutliers (real.stats, 'sackin', signif=0.1^3)  # 11
hist (real.stats$sackin, main='Sackin')
real.stats <- dropOutliers (real.stats, 'colless', signif=0.1^3)  # 11
hist (real.stats$sackin, main='Colless')
# for each gravity stat
psvs <- which (grepl ('psv', colnames (real.stats)))
for (i in psvs) {
  real.stats <- dropOutliers (real.stats, colnames (real.stats)[i],
                              signif=0.1^3)
}
gammas <- which (grepl ('gamma', colnames (real.stats)))
for (i in gammas) {
  real.stats <- dropOutliers (real.stats, colnames (real.stats)[i],
                              signif=0.1^3)
}

# QUICK STATS
# how many tips?
ntrees <- data.frame(n=real.stats$ntaxa[real.stats$ntaxa <= 500 &
                                          !is.na(real.stats$ntaxa)])
mean (ntrees$n)
sd (ntrees$n)
quantile(ntrees$n)
p <- ggplot(ntrees, aes(n)) + geom_histogram() +
  xlab("N. tips") + ylab("N. trees") + theme_bw() +
  theme(text=element_text(size=25))
pdf (file.path (res.dir, 'ntips_histogram.pdf'), 14, 14)
print(p)
dev.off()
# how many polys?
sum (real.stats$poly)
sum (real.stats$poly) *100 /nrow(real.stats)
# how many with bls?
sum (real.stats$bl)
sum (real.stats$bl) *100 /nrow(real.stats)
# when were they published?
mean (real.stats$date)
sd (real.stats$date)
# how many were ultrametric?
sum (real.stats$ultra)
sum (real.stats$ultra) * 100 / nrow (real.stats)
# how were successfully made ultrametric?
psvs <- which (grepl ('psv', colnames (real.stats)))
for (i in psvs) {
  res <-  sum (!is.na (real.stats[ ,i]))
  res.p <- signif (res * 100 / nrow (real.stats), 2)
  name <- colnames (real.stats)[i]
  cat ('For [', name,']: ', res, ' (', res.p, '%)', ']\n', sep='')
}
# how many ..taxonomic..?
length (unique (real.stats$phylum)) - 1
length (unique (real.stats$class)) - 1
length (unique (real.stats$order)) - 1
sum (real.stats$phylum == 'Unknown' & real.stats$class == 'Unknown' &
     real.stats$order == 'Unknown')
unknowns <- which (real.stats$phylum == 'Unknown' & real.stats$class == 'Unknown' &
         real.stats$order == 'Unknown')
real.stats[unknowns, 'title']

# EMPIRICAL
# colless
round (mean (real.stats$colless, na.rm=TRUE), 2)
round (sd (real.stats$colless, na.rm=TRUE), 2)
# sackin
round (mean (real.stats$sackin, na.rm=TRUE), 2)
round (sd (real.stats$sackin, na.rm=TRUE), 2)
# gamma
t.test (real.stats$gamma.pathD8, real.stats$gamma.chronopl)
t.test (real.stats$gamma.pathD8, real.stats$gamma.chronoMPL)
mean (real.stats$gamma.pathD8, na.rm=TRUE)
sd (real.stats$gamma.pathD8, na.rm=TRUE)
mean (real.stats$gamma.chronopl, na.rm=TRUE)
sd (real.stats$gamma.chronopl, na.rm=TRUE)
mean (real.stats$gamma.chronoMPL, na.rm=TRUE)
sd (real.stats$gamma.chronoMPL, na.rm=TRUE)
# still +ve for sourced ultra trees
mean (real.stats$gamma.pathD8[real.stats$ultra], na.rm=TRUE)
mean (real.stats$gamma.chronoMPL[real.stats$ultra], na.rm=TRUE)
mean (real.stats$gamma.chronopl[real.stats$ultra], na.rm=TRUE)
mean (real.stats$psv.pathD8[real.stats$ultra], na.rm=TRUE)
mean (real.stats$psv.chronoMPL[real.stats$ultra], na.rm=TRUE)
mean (real.stats$psv.chronopl[real.stats$ultra], na.rm=TRUE)
# PSV
t.test (real.stats$psv.pathD8, real.stats$psv.chronopl)
t.test (real.stats$psv.pathD8, real.stats$psv.chronoMPL)
mean (real.stats$psv.pathD8, na.rm=TRUE)
sd (real.stats$psv.pathD8, na.rm=TRUE)
mean (real.stats$psv.chronopl, na.rm=TRUE)
sd (real.stats$psv.chronopl, na.rm=TRUE)
mean (real.stats$psv.chronoMPL, na.rm=TRUE)
sd (real.stats$psv.chronoMPL, na.rm=TRUE)

# SIMULATED
table (stats$scenario)
# colless
round (tapply (stats$colless, stats$scenario, mean, na.rm=TRUE), 2)
round (tapply (stats$colless, stats$scenario, sd, na.rm=TRUE) , 2)
round (tapply (extreme$colless, extreme$scenario, mean, na.rm=TRUE), 2)
round (tapply (extreme$colless, extreme$scenario, sd, na.rm=TRUE) , 2)
# sackin
round (tapply (stats$sackin, stats$scenario, mean, na.rm=TRUE), 2)
round (tapply (stats$sackin, stats$scenario, sd, na.rm=TRUE), 2)
round (tapply (extreme$sackin, extreme$scenario, mean, na.rm=TRUE), 2)
round (tapply (extreme$sackin, extreme$scenario, sd, na.rm=TRUE) , 2)
# gamma
res <- tapply (stats$gamma, stats$scenario, mean, na.rm=TRUE)
(res[1]*100/mean (res[2:4]))-100  # % lower DE gamma
t.test (x=stats$gamma[stats$scenario=='DE'],
        y=stats$gamma[stats$scenario!='DE'],
        alternative='less')
round (res, 2)
round (tapply (stats$gamma, stats$scenario, sd, na.rm=TRUE), 2)
real.res <- mean (real.stats$gamma.pathD8, na.rm=TRUE)
(res[1]*100/real.res)-100  # % lower DE gamma
t.test (x=stats$gamma[stats$scenario=='DE'],
        y=real.stats$gamma.pathD8,
        alternative='less')
round (tapply (extreme$gamma, extreme$scenario, mean, na.rm=TRUE), 2)
round (tapply (extreme$gamma, extreme$scenario, sd, na.rm=TRUE) , 2)
round (real.res, 2)
round (sd (real.stats$gamma, na.rm=TRUE), 2)
# PSV
t.test (real.stats$psv.pathD8, real.stats$psv.chronopl)
t.test (real.stats$psv.pathD8, real.stats$psv.chronoMPL)
round (tapply (stats$psv, stats$scenario, mean, na.rm=TRUE), 2)
round (tapply (stats$psv, stats$scenario, sd, na.rm=TRUE), 2)
round (tapply (extreme$psv, extreme$scenario, mean, na.rm=TRUE), 2)
round (tapply (extreme$psv, extreme$scenario, median, na.rm=TRUE), 2)
round (tapply (extreme$psv, extreme$scenario, sd, na.rm=TRUE) , 2)
round (mean (real.stats$psv, na.rm=TRUE), 2)
round (sd (real.stats$psv, na.rm=TRUE), 2)
# age
res <- tapply (stats$age, stats$scenario, mean, na.rm=TRUE)
(res[1]*100/mean (res[2:4]))-100  # % lower DE age
t.test (x=stats$age[stats$scenario=='DE'],
        y=stats$age[stats$scenario!='DE'],
        alternative='less')
round (tapply (extreme$age, extreme$scenario, mean, na.rm=TRUE), 2)
round (tapply (extreme$age, extreme$scenario, sd, na.rm=TRUE) , 2)
round (res, 2)
round (tapply (stats$age, stats$scenario, sd, na.rm=TRUE), 2)


# Compare taxonomic groups
# test for increasing variance between taxonomic groups the lower the taxonomic rank
# first remove unknowns, also only use groups with more than min.sample
# in variance at lower taxonomic levels
taxo.stats <- removeUnknown (real.stats, min.sample=10)
taxo.stats <- removeUnknown (taxo.stats, min.sample=10)  # run twice, taxonomic error
nrow (taxo.stats)
table (taxo.stats$phylum)
table (taxo.stats$class)
table (taxo.stats$order)
# first check the accuracy of the taxonomise script
taxo.stats[sample (1:nrow(taxo.stats), 1), c ('phylum', 'class', 'order')]
# sackin
sackin.phylum <- withinPermTest (taxo.stats$sackin, taxo.stats$phylum)
sackin.class <- withinPermTest (taxo.stats$sackin, taxo.stats$class)
sackin.order <- withinPermTest (taxo.stats$sackin, taxo.stats$order)
sackin.phylum.class <- acrossPermTest (taxo.stats$sackin, taxo.stats$class,
                                       taxo.stats$phylum)
sackin.phylum.order <- acrossPermTest (taxo.stats$sackin, taxo.stats$order,
                                       taxo.stats$phylum)
# colless
colless.phylum <- withinPermTest (taxo.stats$colless, taxo.stats$phylum)
colless.class <- withinPermTest (taxo.stats$colless, taxo.stats$class)
colless.order <- withinPermTest (taxo.stats$colless, taxo.stats$order)
colless.phylum.class <- acrossPermTest (taxo.stats$colless, taxo.stats$class,
                                       taxo.stats$phylum)
colless.phylum.order <- acrossPermTest (taxo.stats$colless, taxo.stats$order,
                                       taxo.stats$phylum)
# psv
psv.pathD8.phylum <- withinPermTest (taxo.stats$psv.pathD8, taxo.stats$phylum)
psv.pathD8.class <- withinPermTest (taxo.stats$psv.pathD8, taxo.stats$class)
psv.pathD8.order <- withinPermTest (taxo.stats$psv.pathD8, taxo.stats$order)
psv.pathD8.phylum.class <- acrossPermTest (taxo.stats$psv.pathD8, taxo.stats$class,
                                        taxo.stats$phylum)
psv.pathD8.phylum.order <- acrossPermTest (taxo.stats$psv.pathD8, taxo.stats$order,
                                        taxo.stats$phylum)
# gamma
gamma.pathD8.phylum <- withinPermTest (taxo.stats$gamma.pathD8, taxo.stats$phylum)
gamma.pathD8.class <- withinPermTest (taxo.stats$gamma.pathD8, taxo.stats$class)
gamma.pathD8.order <- withinPermTest (taxo.stats$gamma.pathD8, taxo.stats$order)
gamma.pathD8.phylum.class <- acrossPermTest (taxo.stats$gamma.pathD8, taxo.stats$class,
                                           taxo.stats$phylum)
gamma.pathD8.phylum.order <- acrossPermTest (taxo.stats$gamma.pathD8, taxo.stats$order,
                                           taxo.stats$phylum)

# FIGURES
# taxonomic
pdf (file.path (res.dir, 'taxonomic.pdf'), 14, 14)
# phyla
instances <- table (real.stats$phylum)
instances <- sort (instances, TRUE)
sum (instances[1:5]) * 100 / sum(instances)  # 85%
real.stats$phylum <- factor (real.stats$phylum, levels = names (instances))
p <- ggplot (real.stats, aes (factor (real.stats$phylum))) +
  geom_bar() + coord_flip() + xlab ('Phylum') + ylab ('N. trees') +
  theme_bw () + theme (text=element_text(size=25))
print (p)
dev.off()
pdf (file.path (res.dir, 'taxonomic_boxplots.pdf'), 20, 14)
p1 <- ggBoxplot (taxo.stats, 'phylum', 'sackin', 'Sackin')
#t.test (x=real.stats$sackin[real.stats$phylum == 'Streptophyta'],
#         y=real.stats$sackin[real.stats$phylum == 'Ascomycota'])
p2 <- ggBoxplot (taxo.stats, 'phylum', 'colless', 'Colless')
p3 <- ggBoxplot (taxo.stats, 'phylum', 'gamma.pathD8', expression(gamma))
p4 <- ggBoxplot (taxo.stats, 'phylum', 'psv.pathD8', 'PSV')
grid.arrange (p1, p2, p3, p4, ncol=2)
# class
instances <- table (taxo.stats$class)
instances <- sort (instances, TRUE)
p1 <- ggBoxplot (taxo.stats, 'class', 'sackin', 'Sackin', 20)
p2 <- ggBoxplot (taxo.stats, 'class', 'colless', 'Colless', 20)
p3 <- ggBoxplot (taxo.stats, 'class', 'psv.pathD8', 'PSV', 20)
p4 <- ggBoxplot (taxo.stats, 'class', 'gamma.pathD8', expression(gamma), 20)
grid.arrange (p1, p2, p3, p4, ncol=2)
# orders
# only use top twenty sampled orders
pull <- taxo.stats$order %in% names (sort (table (taxo.stats$order), TRUE))[1:20]
p1 <- ggBoxplot (taxo.stats[pull, ], 'order', 'sackin', 'Sackin', 20)
p2 <- ggBoxplot (taxo.stats[pull, ], 'order', 'colless', 'Colless', 20)
p3 <- ggBoxplot (taxo.stats[pull, ], 'order', 'psv.pathD8', 'PSV', 20)
p4 <- ggBoxplot (taxo.stats[pull, ], 'order', 'gamma.pathD8', expression(gamma), 20)
grid.arrange (p1, p2, p3, p4, ncol=2)
#t.test (x=real.stats$gamma.pathD8[real.stats$class == 'Actinopterygii'],
#         y=real.stats$gamma.pathD8[real.stats$class == 'Aves'])
# ggBoxplot (real.stats[real.stats$phylum=='Streptophyta', ], 'class', 'psv', 'PSV', 5)
dev.off()
#real.stats[sample (1:nrow(real.stats), 10),c('title', 'phylum')]

# sanity check
stat.names <- c ("colless", "sackin", "psv", "gamma",
                 "age", "pd", 'ntaxa')
pdf(file.path (res.dir, 'sanity_checks.pdf'))
for (stat.name in stat.names) {
  hist (stats[, stat.name], xlab=stat.name, main='Simulated')
  hist (real.stats[, stat.name], xlab=stat.name, main='Empirical')
}
dev.off ()

# Z-scores for simulated trees
pdf (file.path (res.dir, 'tp.pdf'), width=20, height=21.875)
p1 <- tilePlot (stats, stats$colless, legend.title='Colless, Z-score')
p2 <- tilePlot (stats, stats$sackin, legend.title='Sackin, Z-score')
p3 <- tilePlot (stats, stats$gamma, legend.title='Gamma, Z-score')
p4 <- tilePlot (stats, stats$psv, legend.title='PSV, Z-score')
p5 <- tilePlot (stats, stats$age, legend.title='Age, Z-score')
grid.arrange (p1, p2, p3, p4, p5, ncol=2)
dev.off ()

# distances to real trees
pdf (file.path (res.dir, 'tp_dist_to_real.pdf'), width=20, height=17.5)
distances <- abs (stats$colless - mean (real.stats$colless, na.rm=TRUE))
p1 <- tilePlot (stats, distances, legend.title='Colless, Z-score')
distances <- abs (stats$sackin - mean (real.stats$sackin, na.rm=TRUE))
p2 <- tilePlot (stats, distances, legend.title='Sackin, Z-score')
distances <- abs (stats$gamma - mean (real.stats$gamma.pathD8, na.rm=TRUE))
p3 <- tilePlot (stats, distances, legend.title='Gamma, Z-score')
distances <- abs (stats$psv - mean (real.stats$psv.pathD8, na.rm=TRUE))
p4 <- tilePlot (stats, distances, legend.title='PSV, Z-score')
grid.arrange (p1, p2, p3, p4, ncol=2)
dev.off ()

# correlation between sig and imbalance
pdf (file.path (res.dir, 'corr_sig_balance.pdf'))
cor (stats$colless, stats$sig)
p <- ggplot (stats, aes (x=sig, y=colless))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab ('Colless') + xlab (expression (sigma)) +
  theme_bw ()
print (p)
cor (stats$sackin, stats$sig)
p <- ggplot (stats, aes (x=sig, y=sackin))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab ('Sackin') + xlab (expression (sigma)) +
  theme_bw ()
print (p)
dev.off ()

# Correlation between gravity and balance?
# probably due to part of the sample space missing -- low gravity + high imbalance
pull <- !is.na(real.stats$psv.pathD8) & !is.na (real.stats$colless)
cor.test (real.stats$colless[pull], real.stats$psv.pathD8[pull])
# plot (real.stats$colless~real.stats$psv)
# plot (stats$colless~stats$psv)
# real.stats[which (real.stats$psv > 0.95), ]
# cor.test (stats$colless, stats$psv)
pull <- stats$scenario == 'PF'
cor.test (stats$colless[pull], stats$psv[pull])
pull <- stats$scenario == 'DE'
cor.test (stats$colless[pull], stats$psv[pull])
pull <- stats$scenario == 'Pan'
cor.test (stats$colless[pull], stats$psv[pull])
pull <- stats$scenario == 'Eph'
cor.test (stats$colless[pull], stats$psv[pull])


# correlation between eps and gravity for simulations of negative sig
pdf (file.path (res.dir, 'corr_eps_gravity.pdf'))
cor (stats$psv[stats$sig < 0], stats$eps[stats$sig < 0])
p <- ggplot (stats[stats$sig < 0, ], aes (x=eps, y=psv))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab ('PSV') + xlab (expression (epsilon)) +
  theme_bw ()
print (p)
cor (stats$gamma[stats$sig < 0], stats$eps[stats$sig < 0])
p <- ggplot (stats[stats$sig < 0, ], aes (x=eps, y=gamma))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab (expression (gamma)) + xlab (expression (epsilon)) +
  theme_bw ()
print (p)
dev.off ()

# PCA
# 
# filtered <- filter (stats, grain=0.1)
# pca (stats, real.stats, stat.names, 'pca.pdf',
#      ignore.chronos=FALSE)
# pca (filtered, real.stats, stat.names, 'pca_filtered.pdf',
#      ignore.chronos=FALSE)


# Tiles of pca res
stat.names <- c ("colless", "sackin", "psv")
grains <- pca2 (stats, real.stats, stat.names, 'pca_grains.pdf',
                ignore.chronos=FALSE)
real <- grains[nrow (grains), ]
sim <- grains[-nrow (grains), ]
d1 <- abs (real$pc1.mean - sim$pc1.mean)
d1 <- d1/max (d1)
d2 <- abs (real$pc2.mean - sim$pc2.mean)
d2 <- d2/max (d2)
sim$d <- d1 + d2
p <- ggplot (sim, aes (x=eps, y=sig)) + geom_tile (aes (fill=d)) +
  scale_fill_gradient2(low='darkred', mid="red", high="white", name='PC distance') +
  labs (x=expression (epsilon), y=expression (sigma)) +
  theme_bw() + theme(axis.title=element_text(size=8), axis.text=element_text(size=6),
                     legend.text=element_text(size=4), legend.position="left",
                     legend.title=element_text(size=5), legend.key.size=unit(0.25, "cm"),
                     legend.margin=unit(0, "cm"))
# final figure resolution
tiff("~/Desktop/figure.tiff", width=9, height=7, units="cm",
     res=1200)
print (p)
dev.off()

# Looking at PCA of extreme scenarios only
stat.names <- c ("colless", "sackin", "psv")
pca.res <- pca (extreme, real.stats, stat.names, other="quality")
# final figure
p <- PCANoColour (pca.res, file.path (res.dir, 'pca_extreme_1.pdf'))
p <- p + geom_text(hjust = 0, nudge_x = 0.05, vjust=1, nudge_y=-0.05, size=2)
p <- p + theme_bw() + theme(text=element_text(size=4),
                            axis.title=element_text(size=8),
                            axis.text=element_text(size=6),
                            legend.text=element_text(size=6),
                            legend.position="top",
                            legend.text=element_text(size=6),
                            legend.title=element_blank(),
                            legend.key=element_blank(),
                            legend.margin=unit(0, "cm"),
                            legend.key.size=unit(0.25, "cm"))
tiff("~/Desktop/figure.tiff", width=9, height=9, units="cm",
     res=1200)
print(p)
dev.off()

# gamma partitioned (pathD8 results only)
pca.res$x <- pca.res$x[pca.res$x$shape %in% c ('pathD8', 'Sim.'), ]
scenario <- pca.res$x$Scenario
pull <- pca.res$x$gamma > 0 & scenario == 'Emp.'
pca.res$x$Scenario[pull] <- 'Emp. (>0)'
pull <- pca.res$x$gamma < 0 & scenario == 'Emp.'
pca.res$x$Scenario[pull] <- 'Emp. (<0)'
pca.res$x$shape <- NA
plotPCA (pca.res, file.path (res.dir, 'pca_extreme_2.pdf'))
# literature/TB partitioned
pull <- scenario == 'Emp.'
pca.res$x$Scenario[pull] <- paste0('Emp. (', pca.res$x$quality[pull], ')')
# only use empirical trees that are well sampled
cnts <- table(pca.res$x$Scenario)
pull <- pca.res$x$Scenario %in% names(cnts)[cnts > 5]
pca.res$x <- pca.res$x[pull, ]
plotPCA (pca.res, file.path (res.dir, 'pca_extreme_3.pdf'))
