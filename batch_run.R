## 11/07/2014
## D.J. Bennett
## Run 1 and 2 multiple times with multiple args

## Shared parameters
min.time.span <- 5
min.size <- 5
seed.n <- 2
time <- 10
sample <- 0.1
birth <- 1.1
death <- 1
## Unique parameters
i.bias <- c ('none', 'PE', 'FP')

## Run
for (i in 1:3) {
  bias <- i.bias[i]
  source ("1_model.R")
  source ("2_analysis.R")
}