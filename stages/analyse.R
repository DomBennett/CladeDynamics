# 16/07/2014
# D.J. Bennett
# Comparing simulated and empirical trees -- Interactive script

# LIBS
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))
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

# QUICK STATS
# how many tips?
mean (real.stats$ntaxa[real.stats$ntaxa <= 500], na.rm= TRUE)
sd (real.stats$ntaxa[real.stats$ntaxa <= 500], na.rm= TRUE)
# how many polys?
sum (real.stats$poly)
sum (real.stats$poly) *100 /nrow(real.stats)
# how many with bls?
sum (real.stats$bl)
sum (real.stats$bl) *100 /nrow(real.stats)
# when were they published?
mean (real.stats$date)
sd (real.stats$date)
# how many are ultrametric?
sum (real.stats$ul)
sum (real.stats$ul) * 100 / nrow (real.stats)
# how many were made ultrametric?
sum (real.stats$chronos)
sum (real.stats$chronos) * 100 / nrow (real.stats)
# how many suitable for gamma stats?
sum (real.stats$ul) + sum (real.stats$chronos)
# how many ..taxonomic..?
length (unique (real.stats$phylum)) - 1
length (unique (real.stats$class)) - 1
length (unique (real.stats$order)) - 1
sum (real.stats$phylum == 'Unknown' & real.stats$class == 'Unknown' &
       real.stats$order == 'Unknown')
unknowns <- which (real.stats$phylum == 'Unknown' & real.stats$class == 'Unknown' &
         real.stats$order == 'Unknown')
real.stats[unknowns, 'title']
# remove all values where gravity metrics shouldn't have been calclated
real.stats$gamma[!(real.stats$ul | real.stats$chronos)] <- NA
sum (!is.na (real.stats$gamma))  # should be about 469

# PARSE
cat ('\nDropping outliers ....')
# remove branch results where psv is greater than 1 -- this is impossible!
psvs <- which (grepl ('psv', colnames (real.stats)))
for (i in psvs) {
  pull <- real.stats[ ,i] > 1 & !is.na (real.stats[ ,i])
  cat (paste0 (colnames (real.stats)[i], " -- ", sum (pull), '\n'))
  real.stats[pull, i] <- NA
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

# CHECK PSV AND GAMMA BETWEEN METHODS
t.test (real.stats$psv.pathD8, real.stats$psv.chronopl)
t.test (real.stats$psv.pathD8, real.stats$psv.chronoMPL)
t.test (real.stats$gamma.pathD8, real.stats$gamma.chronopl)
t.test (real.stats$gamma.pathD8, real.stats$gamma.chronoMPL)
sd (real.stats$gamma.pathD8, na.rm=TRUE)
sd (real.stats$gamma.chronopl, na.rm=TRUE)
sd (real.stats$gamma.chronoMPL, na.rm=TRUE)

# TABLES
# S2
# colless
round (tapply (stats$colless, stats$scenario, mean, na.rm=TRUE), 2)
round (tapply (stats$colless, stats$scenario, sd, na.rm=TRUE) , 2)
round (tapply (extreme$colless, extreme$scenario, mean, na.rm=TRUE), 2)
round (tapply (extreme$colless, extreme$scenario, sd, na.rm=TRUE) , 2)
round (mean (real.stats$colless, na.rm=TRUE), 2)
round (sd (real.stats$colless, na.rm=TRUE), 2)
# sackin
round (tapply (stats$sackin, stats$scenario, mean, na.rm=TRUE), 2)
round (tapply (stats$sackin, stats$scenario, sd, na.rm=TRUE), 2)
round (tapply (extreme$sackin, extreme$scenario, mean, na.rm=TRUE), 2)
round (tapply (extreme$sackin, extreme$scenario, sd, na.rm=TRUE) , 2)
round (mean (real.stats$sackin, na.rm=TRUE), 2)
round (sd (real.stats$sackin, na.rm=TRUE), 2)
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
round (tapply (stats$psv, stats$scenario, mean, na.rm=TRUE), 2)
round (tapply (stats$psv, stats$scenario, sd, na.rm=TRUE), 2)
round (tapply (extreme$psv, extreme$scenario, mean, na.rm=TRUE), 2)
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
# increasing variance between taxonomic groups the lower the taxonomic rank
sackin.phylum <- tapply (real.stats$sackin, factor (real.stats$phylum), mean, na.rm=TRUE)
var (sackin.phylum)
sackin.class <- tapply (real.stats$sackin, factor (real.stats$class), mean, na.rm=TRUE)
var (sackin.class)
sackin.order <- tapply (real.stats$sackin, factor (real.stats$order), mean, na.rm=TRUE)
var (sackin.order)
psv.phylum <- tapply (real.stats$psv, factor (real.stats$phylum), mean, na.rm=TRUE)
var(psv.phylum, na.rm=TRUE)
psv.class <- tapply (real.stats$psv, factor (real.stats$class), mean, na.rm=TRUE)
var(psv.class, na.rm=TRUE)
psv.order <- tapply (real.stats$psv, factor (real.stats$order), mean, na.rm=TRUE)
var(psv.order, na.rm=TRUE)
# increase for both balance and gravity
var (sackin.order)*100/var (sackin.phylum)  # 232% increase in Sackin
var (psv.order, na.rm=TRUE)*100/var (psv.phylum, na.rm=TRUE)  # 127% increase in Sackin

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
ggBoxplot (real.stats, 'phylum', 'sackin', 'Sackin')
# t.test (x=real.stats$sackin[real.stats$phylum == 'Streptophyta'],
#         y=real.stats$sackin[real.stats$phylum == 'Ascomycota'])
ggBoxplot (real.stats, 'phylum', 'colless', 'Colless')
ggBoxplot (real.stats, 'phylum', 'gamma', expression(gamma))
ggBoxplot (real.stats, 'phylum', 'psv', 'PSV')
# class
instances <- table (real.stats$class)
instances <- sort (instances, TRUE)
ggBoxplot (real.stats, 'class', 'sackin', 'Sackin', 20)
ggBoxplot (real.stats, 'class', 'colless', 'Colless', 20)
ggBoxplot (real.stats, 'class', 'psv', 'PSV', 20)
ggBoxplot (real.stats, 'class', 'gamma', expression(gamma), 20)
# orders
ggBoxplot (real.stats, 'order', 'sackin', 'Sackin', 20)
ggBoxplot (real.stats, 'order', 'colless', 'Colless', 20)
ggBoxplot (real.stats, 'order', 'psv', 'PSV', 20)
ggBoxplot (real.stats, 'order', 'gamma', expression(gamma), 20)
# t.test (x=real.stats$gamma[real.stats$class == 'Actinopterygii'],
#         y=real.stats$gamma[real.stats$class == 'Aves'])
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
pdf (file.path (res.dir, 'tp.pdf'), width=8)
p <- tilePlot (stats, stats$colless, legend.title='Colless, Z-score')
print (p)
p <- tilePlot (stats, stats$sackin, legend.title='Sackin, Z-score')
print (p)
p <- tilePlot (stats, stats$gamma, legend.title='Gamma, Z-score')
print (p)
p <- tilePlot (stats, stats$psv, legend.title='PSV, Z-score')
print (p)
p <- tilePlot (stats, stats$age, legend.title='Age, Z-score')
print (p)
dev.off ()

# distances to real trees
pdf (file.path (res.dir, 'tp_dist_to_real.pdf'), width=8)
distances <- abs (stats$colless - mean (real.stats$colless, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='Colless, Z-score')
print (p)
distances <- abs (stats$sackin - mean (real.stats$sackin, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='Sackin, Z-score')
print (p)
distances <- abs (stats$gamma - mean (real.stats$gamma, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='Gamma, Z-score')
print (p)
distances <- abs (stats$psv - mean (real.stats$psv.chronopl, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='PSV, Z-score')
print (p)
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
pull <- !is.na(real.stats$psv) & !is.na (real.stats$colless)
cor.test (real.stats$colless[pull], real.stats$psv[pull])
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
stat.names <- c ("colless", "sackin", "psv")
filtered <- filter (stats, grain=0.1)
pca (stats, real.stats, stat.names, 'pca.pdf',
     ignore.chronos=FALSE)
pca (filtered, real.stats, stat.names, 'pca_filtered.pdf',
     ignore.chronos=FALSE)
grains <- pca2 (stats, real.stats, stat.names, 'pca_grains.pdf',
                ignore.chronos=FALSE)

# Tiles of pca res
pdf (file.path (res.dir, 'tp_pca.pdf'))
real <- grains[nrow (grains), ]
sim <- grains[-nrow (grains), ]
d1 <- abs (real$pc1.mean - sim$pc1.mean)
d1 <- d1/max (d1)
d2 <- abs (real$pc2.mean - sim$pc2.mean)
d2 <- d2/max (d2)
sim$d <- d1 + d2
p <- ggplot (sim, aes (x=eps, y=sig)) + geom_tile (aes (fill=d)) +
  scale_fill_gradient2(low='red', mid='red', high='white', name='PC distance') +
  labs (x=expression (epsilon), y=expression (sigma)) +
  theme_bw() + theme (axis.title=element_text(size=25))
print (p)
dev.off()

# Looking at PCA of extreme scenarios only
stat.names <- c ("colless", "sackin", "psv")
pca.res <- pca (extreme, real.stats, stat.names)
plotPCA (pca.res, file.path (res.dir, 'pca_extreme_1.pdf'))
# gamma partitioned (pathD8 results only)
pca.res$x <- pca.res$x[pca.res$x$shape %in% c ('pathD8', 'Sim.'), ]
emp.gamma <- pca.res$x$gamma > 0 & pca.res$x$shape == 'pathD8'
pca.res$x$shape <- ifelse (emp.gamma, 'High', 'Low')
plotPCA (pca.res, file.path (res.dir, 'pca_extreme_2.pdf'))
