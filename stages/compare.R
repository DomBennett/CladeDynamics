## 16/07/2014
## D.J. Bennett
## Comparing simulated and empirical trees

## Libraries
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))

## Dirs
data.dir <- file.path ('data', 'treestats')

## Parameters
analysis.name <- 'analysis_4'
data.dir <- file.path ('data', 'treestats')
empirical.file <- 'min50_max500.Rd'

## Input and generation
res.dir <- file.path ('results', analysis.name)
runlog <- file.path (res.dir, 'runlog.csv')
# get metadata
metadata <- read.csv (runlog, stringsAsFactors=FALSE)
if (!file.exists (file.path (res.dir, 'stats.Rd'))) {
  # get simulated trees and calc stats
  trees <- readTrees (metadata, res.dir, runlog)
  stats <- calcTreeStats(trees)
  stats <- cbind (metadata, stats)
  stats <- getScenarios(stats)
  ed.values <- getEDs (trees, stats$scenario)
  save (stats, ed.values, file=file.path (res.dir, 'stats.Rd'))
} else {
  load (file.path (res.dir, 'stats.Rd'))
}
# load pre-calculated empirical tree stats -- real.stats and real.ed.values
load (file.path (data.dir, empirical.file))

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
tapply (stats$gamma, stats$scenario, mean, na.rm=TRUE)
tapply (stats$gamma, stats$scenario, sd, na.rm=TRUE)
mean (real.stats$gamma, na.rm=TRUE)
sd (real.stats$gamma, na.rm=TRUE)
# PSV
tapply (stats$psv, stats$scenario, mean, na.rm=TRUE)
tapply (stats$psv, stats$scenario, sd, na.rm=TRUE)
tapply (real.stats$psv, real.stats$ul, mean, na.rm=TRUE)
tapply (real.stats$psv, real.stats$ul, sd, na.rm=TRUE)
# age
tapply (stats$age, stats$scenario, mean, na.rm=TRUE)
tapply (stats$age, stats$scenario, sd, na.rm=TRUE)
tapply (real.stats$age, real.stats$ul, mean, na.rm=TRUE)
tapply (real.stats$age, real.stats$ul, sd, na.rm=TRUE)

## Figures
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
p <- ggplot (stats, aes (x=sig, y=colless))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab ('Colless') + xlab (expression (sigma)) +
  theme_bw ()
print (p)
p <- ggplot (stats, aes (x=sig, y=sackin))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab ('Sackin') + xlab (expression (sigma)) +
  theme_bw ()
print (p)
dev.off ()

# figure 6 -- correlation between eps and loading for simulations of negative sig
pdf (file.path (res.dir, 'figure_6.pdf'))
p <- ggplot (stats[stats$sig < 0, ], aes (x=eps, y=psv))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab ('PSV') + xlab (expression (epsilon)) +
  theme_bw ()
print (p)
p <- ggplot (stats[stats$sig < 0, ], aes (x=eps, y=gamma))
p <- p + geom_point () + stat_smooth (method='lm') +
  ylab (expression (gamma)) + xlab (expression (epsilon)) +
  theme_bw ()
print (p)
dev.off ()

# figure 7 -- ED distributions
pdf (file.path (res.dir, 'figure_7.pdf'))
#ed.data <- rbind (ed.values, real.ed.values)  # TODO
ed.data <- ed.values
# TODO -- plot each of the scenarios with real EDs
p <- ggplot (ed.data, aes (x=ed.values, fill=groups))
p <- p + geom_density (alpha=0.5) + xlab ('ED, Z-score') + theme_bw()
print (p)
dev.off ()

# figure 8 -- PCA
stat.names <- c ("colless", "sackin", "psv")
pca (stats, real.stats, stat.names, 'figure8_withoutchronos.pdf',
     ignore.chronos=TRUE)
pca (stats, real.stats, stat.names, 'figure8_withchronos.pdf',
     ignore.chronos=FALSE)
