# 19/01/2015
# How does gamma change with respect to psi?

# LIBS
source (file.path ('tools', 'model_tools.R'))

# PARAMETERS
n <- 100


# INIT
trees <- list ()
psis <- runif (n, -1, 1)
gammas <- sizes <- ages <- bts <- rep (NA, n)
# use a random tree of 100 taxa
seed.tree <- compute.brlen (rtree (100))
# need to change labels to ensure no overlap with tip labels in runEDBMM()
seed.tree$tip.label <- paste0 ('t', 1:getSize (seed.tree))
seed.tree$edge.length <- seed.tree$edge.length*0.1
getSize (seed.tree, 'rtt')
gammaStat (seed.tree)

# LOOP
for (i in 1:n) {
  cat ('\n.... [', i, '/', n, ']', sep = '')
  # model with equal birth and death to remove its impact
  tree <- runEDBMM (birth = 1, death = 1, stop.at = 2,
                    stop.by = 't', psi = psis[i],
                    seed.tree = seed.tree)
  trees <- c (trees, list (tree))
  gammas[i] <- gammaStat (tree)
  sizes[i] <- getSize (tree)
  ages[i] <- getSize (tree, 'rtt')
  bts[i] <- mean (branching.times(tree))
}

# PLOTING
#plot (gammas ~ bds)
#plot (sizes ~ bds)
#plot (ages ~ bds)
plot (gammas ~ psis)  # all gammas are negative, must be a hangover from the rtree
plot (bts ~ psis)
plot (sizes ~ psis)  # no relationship for these
plot (ages ~ psis)

# HIGH PSI --> LOW GAMMA
# LOW PSI --> HIGH GAMMMA
# This makes sense as negative psi means recently speciated are more likely to speciate,
# concentrating nodes towards the present.

# REPEAT WO RTREE
n <- 100
trees <- list ()
psis <- runif (n, -1, 1)
gammas <- sizes <- ages <- bts <- rep (NA, n)

# LOOP
for (i in 1:n) {
  cat ('\n.... [', i, '/', n, ']', sep = '')
  tree <- runEDBMM (birth = 20, death = 10, stop.at = 100,
                    stop.by = 'n', psi = psis[i])
  trees <- c (trees, list (tree))
  gammas[i] <- gammaStat (tree)
  sizes[i] <- getSize (tree)
  ages[i] <- getSize (tree, 'rtt')
  bts[i] <- mean (branching.times(tree))
}

# PLOTTING
plot (gammas ~ psis)
plot ((bts/ages) ~ psis)
plot (bts ~ ages)
plot (sizes ~ psis)  # no relationship for these
plot (ages ~ psis)