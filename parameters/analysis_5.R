# D.J. Bennett
# 10,000 simulations

## Parameter descriptions
# n.model -- the number of trees to simulate
# seed -- starting number of taxa in random seed tree (must be >= 2)
# stop.by -- aim for number of taxa (n) or amount of time (t) in end tree
# max.ntaxa -- max ntaxa in a tree
# min.ntaxa -- min ntaxa in a tree
# birth -- how many births per unit of branch length?
# death -- how many deaths per unit of branch length?
# bias -- what type of ED? 'PE', 'ES' or 'FP'
# min.sig -- min power determing the effect on speciation
# max.sig -- max power determing the effect on speciation
# min.eps -- min power determing the effect on extinction
# max.eps -- max power determing the effect on extinction
# reference -- normalise shape stats with a Yule reference?
# iterations -- n trees in Yule distribution if reference is True

## Parameter declarations (list of lists)
analysis.5 <- list (n.model = 10000, seed = 2,
                    max.birth = 2, min.birth = 2,
                    max.death = 1, min.death = 1,
                    bias = 'FP', stop.by = 'n',
                    max.ntaxa = 500, min.ntaxa = 50,
                    min.sig = -1, max.sig = 1,
                    min.eps = -1, max.eps = 1,
                    reference = TRUE,
                    iterations = 100)
analysis.parameters <- list (analysis_5=analysis.5)