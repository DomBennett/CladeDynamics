# Evolutionary Distinctiveness Biased Markov Modelling (EDBMM)
A R pipeline for exploring how an evolutionary distinctiveness bias
in a tree-growth Markov model affects tree shape and clade dynamics.

## System requirements
*R version
  +3.0.1+
*Deps
  +`plyr`
  +`ape`
  +`geiger`
  +`apTreeshape`
  +`ggplot2`
  +`MoreTreeTools`
*Optional
  +`test_that`

##Pipeline structure
###Scripts
*`batch_run.R`: run all analyses
*`run_analysis_1.R`: run analysis 1
*`run_analysis_2.R`: run analysis 2
*`run_analysis_3.R`: run analysis 3
*`run_analysis_4.R`: run analysis 4
*`run_tests.R`: run this script to test all custom functions using the test_that package
###Folders
*data: contains real phylogenetic trees and their stats, see its README
*original_results: the results as used in publication
*tools: contains custom functions associated with each script
*other: contains scripts for running pre-calculations, see its README
*stages: scripts for running differents steps of the analysis pipeline
*sanity_checks: scripts for testing assumptions and ideas are good
*results: this folder will hold all the results generated for each analysis
  +analysis_1
  +analysis_2
  +analysis_3
  +analysis_4

##Notes
1.`MoreTreeTools` is not available from CRAN but can be installed with
the package `devtools` via GitHub:
```{r}
library (devtools)
install_github ('DomBennett/MoreTreeTools')
```
2. Results can be reproduce results it in an R session:
```{r}
source ('batch_run.R')
```
Or, via command-line:
```
Rscript batch_run.R
```
3. The original results took ## hours and ## minutes to produce running on ...

##Author
D.J.Bennett (ICL & ZSL)