## 18/07/2014
## D.J. Bennett
R version:
 R 3.0.1+
Deps:
 plyr, ape, geiger, apTreeshape, ggplot2, MoreTreeTools
Optional:
 test_that
Scripts:
 '1_model.R': model tree growth
 '2_analyse.R': produce report graphs and stats for output from '1_model.R'
 'run.R': run '1_model.R' and '2_analyse.R' in succession
 '3_compare.R': compare stats from multiple model runs and real trees
 'batch_run.R': run all scripts with parameters specified in 'parameters.R' to model 
  different tree growth stats
 'run_test.R': run this script to test all custom functions using the test_that package
Folders:
 data: contains real phylogenetic trees, see its README
 original_results: the results as used in publication
 tools: contains custom functions associated with each script
Notes:
 1. MoreTreeTools is not available from CRAN but can be installed with the package
 devtools and via GitHub:
  > library (devtools)
  > install_github ('DomBennett/MoreTreeTools')
 2. Results can be reproduce results it in an R session:
   > source ('batch_run.R')
  Or, via command-line:
   > Rscript batch_run.R
 3. The original results took ## hours and ## minutes to produce running on 
