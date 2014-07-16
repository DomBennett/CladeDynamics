## 16/07/2014
## D.J. Bennett
## Analysis of different modelled rises and falls of clades

## Libraries
source (file.path ('tools', 'analysis_tools.R'))

## Dirs
# read in last run folder
res.dir <- read.delim (file.path ('results', 'run_log.txt'),
                       header = FALSE, stringsAsFactors = FALSE)[ ,1]
# read in last n.compare runs
res.dirs <- res.dir[(length (res.dir) - n.compare):length (res.dir)]

## Read in results + plot ages
filename <- file.path (
  'results', paste0 (
    'cladeages_for_',paste(res.dirs, collapse = '_'), '.pdf'))
pdf (filename)
histClageAges (res.dirs)
dev.off ()