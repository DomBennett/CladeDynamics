## 07/08/2014
## D.J. Bennett
## Model trees using EDBMM

## Libraries
source (file.path ('tools', 'model_tools.R'))

## Dirs
# create a unique dir based on time for every run
counter <- nrow (read.csv (runlog))
# create a tree filename
treefilename <- paste0 ('tree', counter, '.tre')
# record parameters in runlog
parameters <- data.frame (treefilename, strength, bias,
                          birth, death)
write.table (parameters, runlog, sep = ',', append = TRUE,
             col.names = FALSE, row.names = FALSE)

## Model
cat (paste0 ('\nModelling tree of size [', stop.at, '] [',
             stop.by, '] ...'))
tree <- runEDBMM (birth = birth, death = death,
                  stop.at = stop.at, stop.by = stop.by,
                  strength = strength, bias = bias,
                  fossils = FALSE, record = record)
## TODO
# get this to work for an analysis where birth and death are equal

## Saving results
cat ('\nSaving results ...')
write.tree (tree, file = file.path (res.dir, treefilename))