## 11/07/2014
## D.J. Bennett
## Run 1 and 2 multiple times with multiple args

## Shared parameters
burnin <- 1 # time to throw away to remove bias of start
min.time.span <- 5
min.size <- 5
seed.n <- 1000
time <- 100
sample <- 0.1
birth <- 1
death <- 1
## Unique parameters
i.bias <- c ('none', 'PE', 'FP', 'iPE', 'iFP')

## Run
for (i in 1:length (i.bias)) {
  bias <- i.bias[i]
  cat (paste0 ('Running for bias: [', bias, ']\n\n'))
  source ("1_model.R", echo = TRUE)
  source ("2_analysis.R", echo = TRUE)
}