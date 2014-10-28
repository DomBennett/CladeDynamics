## D.J. Bennett
## 21/10/2014
## What do clade rises and falls looks like?

## Sources
source (file.path ('tools', 'model_tools.R'))
source (file.path ('tools', 'clade_tools.R'))

## Model trees
t.stop <- 20
sample <- t.stop*.01
seed.tree <- compute.brlen (rtree (100))
seed.tree$tip.label <- paste0 ('t', 1:getSize (seed.tree))
seed.tree$node.label <- paste0 ('n', 1:seed.tree$Nnode)
trees <- runEDBMM (birth = 1, death = 1, stop.at = t.stop, stop.by = 't',
                   bias = 'FP',record = TRUE, fossils = FALSE, psi = -1,
                   progress.bar = 'time', sample = sample, seed.tree = seed.tree)
clades <- getCladeSuccess (trees, 1)
clade.stats <- calcCladeStats (clades)
# what's the mean time span, cm and cg for clades
# -- over 50 spp
# -- appeared and disappeared during simulation
bool <- clade.stats$max.size > 50 & clade.stats$start != 1 & clade.stats$end != length (trees)
sum (bool)
clade.stats[bool, ]
mean(clade.stats[bool,]$cm)

## Exploring speeding up runEDBMM
Rprof("out.out")
t.stop <- 10
sample <- t.stop*.1
trees <- runEDBMM (birth = 2, death = 1, stop.at = 20, stop.by = 't',
                   bias = 'FP',record = TRUE, fossils = FALSE, psi = -1,
                   progress.bar = 'none', sample = sample)
Rprof(NULL)
summaryRprof("out.out")
