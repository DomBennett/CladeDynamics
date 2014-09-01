## 01/09/2014
## D.J. Bennett
## Exploring what the natural tree stats look like

load ('data/natural_tree_stats_t100_l10.Rd')

# bootstrap 95% confidence intervals
data <- natural.tree.stats$colless.stat
reps <- 999
sampled.means <- rep (NA, reps)
for (i in 1:reps) {
  sampled.means[i] <- mean (sample (data, 100, replace = TRUE))
}
hist (sampled.means)
mean (sampled.means)
quantile (sampled.means, 0.95)
quantile (sampled.means, 0.05)
