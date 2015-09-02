# Evolutionary Distinctness Biased Markov Model (EDBMM)

An R pipeline for exploring how an evolutionary distinctness bias in a
tree-growth Markov model affects tree shape and clade dynamics in order to test
the reality of the living fossil.

Data and results files are not provided in this repository, only the code is.
The final set of files and folders used for publication can be found [here]().
 (*Publication pending*)

**System requirements**

Run the `install_deps.R` script to install all dependent packages automatically.

* OS
  + UNIX (but readily modifiable for Windows)
* R version
  + 3+
* R packages
  + `plyr`
  + `ape`
  + `geiger`
  + `apTreeshape`
  + `ggplot2`
  + `MoreTreeTools`
  + `outliers`
  + `doMC`
  + `foreach`
  + `test_that` (optional)

**Directory structure**

```
-- data/
---- raw_trees/
------ literature/
-------- [manually added]
------ treebase/
---- parsed_trees/
---- treestats/
-- stages/
---- [all stage .R scripts]
-- tools/
---- [all tool .R scripts]
-- results/
----- [any folders named by analysis as specified in run.R]
-- other/
-- sanity_checks/
```

**Pipeline**

The pipeline works by calling either `setup.R` or `run.R`. These scripts call
stage scripts which can be found in [`/stages`](https://github.com/DomBennett/Project-EDBMM/tree/master/stages).
The stages scripts depend on custom functions found in the [`/tools`](https://github.com/DomBennett/Project-EDBMM/tree/master/tools) folder.

**setup.R**

This phase of the pipeline sources and calculates statistics from empirical
trees:

1. Download trees from TreeBase
2. Parse trees
3. Calculate tree shape statistics

**run.R**

This reproduces all the results:

1. Model trees according to parameters in `run.R`
2. Compare results from modelled trees with statistics generated at setup.

N.B. the taxonomise and clade stages must be run separately, these are not part
of the pipeline and are were post-hoc.

**Testing**

Run `test.R` to make sure core functions are working.

**Author**

Dom Bennett
