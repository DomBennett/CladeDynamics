## 07/08/2014
## D.J. Bennett
## Model trees for set number of taxa


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
cat (paste0 ('\nModelling tree of size [', n.taxa, '] ...'))
tree <- growMRMMTree (birth = birth, death = death,
                      stop.at = n.taxa, stop.by = 'max.n',
                      strength = strength, bias = bias,
                      fossils = FALSE)

## Saving results
cat ('\nSaving results ...')
write.tree (tree, file = file.path (res.dir, treefilename))