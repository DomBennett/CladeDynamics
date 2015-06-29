## 16/07/2014
## D.J. Bennett
## Comparing simulated and empirical trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))
library (outliers)

## Dirs
data.dir <- file.path ('data', 'treestats')

## Parameters
analysis.name <- 'analysis_5'
res.dir <- file.path ('results', analysis.name)
data.dir <- file.path ('data', 'treestats')
empirical.file <- 'min50_max500.Rd'

## Input and generation
stats <- readIn (analysis.name)
# load pre-calculated empirical tree stats -- real.stats and real.ed.values
load (file.path (data.dir, empirical.file))
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


# Quick stats
# how many polys?
sum (real.stats$poly)
sum (real.stats$poly) *100 /nrow(real.stats)
# how many with bls?
sum (real.stats$bl)
sum (real.stats$bl) *100 /nrow(real.stats)
# when were they published?
mean (real.stats$date)
sd (real.stats$date)
# how many tips?
mean (real.stats$ntaxa[real.stats$ntaxa <= 500], na.rm= TRUE)
sd (real.stats$ntaxa[real.stats$ntaxa <= 500], na.rm= TRUE)
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

# remove all values where gravity metrics shouldn't have been calclated
real.stats$gamma[!(real.stats$ul | real.stats$chronos)] <- NA
sum (!is.na (real.stats$gamma))  # should be about 490

# remove branch results where psv is greater than 1 -- this is impossible!
pull <- real.stats$psv > 1
sum (pull)  # 34 lost
real.stats$psv[pull] <- NA
real.stats$gamma[pull] <- NA

# Extremely conservative removal of outliers
cat ('\nDropping outliers ....')
real.stats <- dropOutliers (real.stats, 'sackin', signif=0.1^3)  # 11
hist (real.stats$sackin, main='Sackin')
real.stats <- dropOutliers (real.stats, 'colless', signif=0.1^3)  # 11
hist (real.stats$sackin, main='Colless')
real.stats <- dropOutliers (real.stats, 'gamma', signif=0.1^3)  # 0
hist (real.stats$gamma, main='Gamma')
real.stats <- dropOutliers (real.stats, 'psv', signif=0.1^3)  # 0
hist (real.stats$psv, main='PSV')

# Table 3
# colless
tapply (stats$colless, stats$scenario, mean, na.rm=TRUE)
tapply (stats$colless, stats$scenario, sd, na.rm=TRUE)
mean (real.stats$colless, na.rm=TRUE)
sd (real.stats$colless, na.rm=TRUE)
# sackin
tapply (stats$sackin, stats$scenario, mean, na.rm=TRUE)
tapply (stats$sackin, stats$scenario, sd, na.rm=TRUE)
mean (real.stats$sackin, na.rm=TRUE)
sd (real.stats$sackin, na.rm=TRUE)
# gamma
res <- tapply (stats$gamma, stats$scenario, mean, na.rm=TRUE)
(res[1]*100/mean (res[2:4]))-100  # % lower DE gamma
t.test (x=stats$gamma[stats$scenario=='DE'],
        y=stats$gamma[stats$scenario!='DE'],
        alternative='less')
tapply (stats$gamma, stats$scenario, sd, na.rm=TRUE)
tapply (real.stats$gamma, real.stats$ul | real.stats$chronos,
        mean, na.rm=TRUE)
tapply (real.stats$gamma, real.stats$ul | real.stats$chronos,
        sd, na.rm=TRUE)
real.res <- mean (real.stats$gamma, na.rm=TRUE)
(res[1]*100/real.res)-100  # % lower DE gamma
t.test (x=stats$gamma[stats$scenario=='DE'],
        y=real.stats$gamma,
        alternative='less')
sd (real.stats$gamma, na.rm=TRUE)
# PSV
tapply (stats$psv, stats$scenario, mean, na.rm=TRUE)
tapply (stats$psv, stats$scenario, sd, na.rm=TRUE)
tapply (real.stats$psv, real.stats$ul, mean, na.rm=TRUE)
tapply (real.stats$psv, real.stats$ul, sd, na.rm=TRUE)
mean (real.stats$psv, na.rm=TRUE)
sd (real.stats$psv, na.rm=TRUE)
# age
res <- tapply (stats$age, stats$scenario, mean, na.rm=TRUE)
(res[1]*100/mean (res[2:4]))-100  # % lower DE age
t.test (x=stats$age[stats$scenario=='DE'],
        y=stats$age[stats$scenario!='DE'],
        alternative='less')
tapply (stats$age, stats$scenario, sd, na.rm=TRUE)
tapply (real.stats$age, real.stats$ul, mean, na.rm=TRUE)
tapply (real.stats$age, real.stats$ul, sd, na.rm=TRUE)

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

## Figures
# taxonomic
pdf (file.path (res.dir, 'taxonomic.pdf'), 14, 14)
# phyla
instances <- table (real.stats$phylum)
instances <- sort (instances, TRUE)
sum (instances[1:5]) * 100 / sum(instances)  # 85%
real.stats$phylum <- factor (real.stats$phylum, levels = names (instances))
p <- ggplot (real.stats, aes (factor (real.stats$phylum))) +
  geom_bar() + coord_flip() + xlab ('Phylum') + ylab ('N. trees') +
  theme (text=element_text(size=25))
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

# figure 3 -- Z-scores for simulated trees
pdf (file.path (res.dir, 'figure_3.pdf'), width=8)
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

# figure 4 -- distances to real trees
pdf (file.path (res.dir, 'figure_4.pdf'), width=8)
distances <- abs (stats$colless - mean (real.stats$colless, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='Colless, Z-score')
print (p)
distances <- abs (stats$sackin - mean (real.stats$sackin, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='Sackin, Z-score')
print (p)
distances <- abs (stats$gamma - mean (real.stats$gamma, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='Gamma, Z-score')
print (p)
distances <- abs (stats$psv - mean (real.stats$psv, na.rm=TRUE))
p <- tilePlot (stats, distances, legend.title='PSV, Z-score')
print (p)
dev.off ()

# figure 5 -- correlation between sig and imbalance
pdf (file.path (res.dir, 'figure_5.pdf'))
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


# figure 6 -- correlation between eps and gravity for simulations of negative sig
pdf (file.path (res.dir, 'figure_6.pdf'))
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

# figure 7 -- ED distributions
# pdf (file.path (res.dir, 'figure_7.pdf'))
# #ed.data <- rbind (ed.values, real.ed.values)  # TODO
# ed.data <- ed.values
# # TODO -- plot each of the scenarios with real EDs
# p <- ggplot (ed.data, aes (x=ed.values, fill=groups))
# p <- p + geom_density (alpha=0.5) + xlab ('ED, Z-score') + theme_bw()
# print (p)
# dev.off ()

# figure 8 -- PCA
stat.names <- c ("colless", "sackin", "psv")
filtered <- filter (stats, grain=0.1)
pca (stats, real.stats, stat.names, 'figure8.pdf',
     ignore.chronos=FALSE)
pca (filtered, real.stats, stat.names, 'figure8_filtered.pdf',
     ignore.chronos=FALSE)
grains <- pca2 (stats, real.stats, stat.names, 'figure8_grains.pdf',
                ignore.chronos=FALSE)

# figure 10 -- tiles of pca res
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

# figure 9 -- looking at PCA of extreme scenarios only
stat.names <- c ("colless", "sackin", "psv")
extreme <- rbind (readIn ('Pan'), readIn ('Eph'),
                  readIn ('DE'), readIn ('PF'))
pca (extreme, real.stats, stat.names, 'figure9.pdf',
     ignore.chronos=FALSE)

# Table 3. test what proportion overlaps with real


## Clade analysis
# add cg and cm stats to stats
min.size <- 50
stats$cg.mean <- stats$cm.mean <- stats$cg.sd <- stats$cm.sd <- NA
clades.files <- sub ('\\.tre', '_clades.csv', stats$treefilename)
cladestats.files <- sub ('\\.tre', '_clade_stats.csv', stats$treefilename)
for (i in 1:nrow (stats)) {
  if (!file.exists (file.path ('results', name,
                               clades.files[i]))) {
    next
  }
  clades <- read.csv (file.path ('results', name,
                                 clades.files[i]))[,-1]
  clade.stats <- read.csv (file.path ('results', name,
                                      cladestats.files[i]))[,-1]
  # filtering...
  # ... ignore clades that started
  clade.stats <- clade.stats[clade.stats$start != 1,]
  # ... ignore clades that were still extant
  clade.stats <- clade.stats[clade.stats$end != max (clade.stats$end),]
  # ... only take clades over a certain size
  clade.stats <- clade.stats[clade.stats$max.size > min.size, ]
  # add to stats
  stats$cg.mean[i] <- mean (clade.stats$cg, na.rm=TRUE)
  stats$cg.sd[i] <- sd (clade.stats$cg, na.rm=TRUE)
  stats$cm.mean[i] <- mean (clade.stats$cm, na.rm=TRUE)
  stats$cm.sd[i] <- sd (clade.stats$cm, na.rm=TRUE)
}
# cg and cm by scenario
tapply (stats$cm.mean, stats$scenario, mean, na.rm=TRUE)
tapply (stats$cg.mean, stats$scenario, mean, na.rm=TRUE)