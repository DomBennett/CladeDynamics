## 05/07/2014
## D.J. Bennett
## Run 1 and 2

## Model Parameters
seed.n <- 2 # how big should the initial tree be?
time <- 5 # how many units of branch length should the tree grow by?
sample <- 0.01 # how often should sampling the success of a tree occur?
birth <- 1.1 # how many births to deaths per unit of branch length?
death <- 1
bias <- 'none' # 'none', 'PE' or 'FP'
## Analysis Parameters
min.time.span = 5 # minimum amount of time a clade exists for it to be normalised plotted
min.size = 5 # minimum maximum size for normalised plotting

## Run
source ("1_model.R")
source ("2_analysis.R")