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
* R packages (version used)
  + `plyr` (1.8.1)
  + `ape` (3.2)
  + `geiger` (2.0.3)
  + `apTreeshape` (1.4.5)
  + `ggplot2` (1.0.0)
  + `MoreTreeTools` (0.0.1)
  + `outliers` (0.14)
  + `doMC` (1.3.3, not for Windows)
  + `foreach` (1.4.2)
  + `test_that` (0.9.1, optional)

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
2. Calculate tree shape statistics

In the this script you will need to specify the parameters in a list structure.
For example the 10,000 trees analysis used in the publication was called
analysis_5 and was simulated with this list object:

```{R}
analysis.5 <- list (n.model = 10000, seed = 2,
                    max.birth = 2, min.birth = 2,
                    max.death = 1, min.death = 1,
                    bias = 'FP', stop.by = 'n',
                    max.ntaxa = 500, min.ntaxa = 50,
                    min.sig = -1, max.sig = 1,
                    min.eps = -1, max.eps = 1,
                    reference = TRUE,
                    iterations = 100)
```

The results from this analysis were subsequently supplemented with
results generated outside of the pipeline structure from an additional set of
analyses looking at the extremes of the scenario parameter space (analysis
  parameters in latest version of run.R) and with the results from the
taxonomise and clade stages [first 1000 results]. All results described in the
publication were then generated with the `compare.R` script.

**Testing**

Run `test.R` to make sure core functions are working.

**Author**

Dom Bennett
