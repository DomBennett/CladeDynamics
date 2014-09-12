## 12/09/2014
## D.J. Bennett
## Tools for downloading trees

## Functions
safeConnect <- function (expr, wait = 3600, trys = 10,
                         verbose = TRUE) {
  ## Safely connect to online datasources using expr
  ##  in large data pipelines without pipeline failing
  # I use an hour wait and ten trys because TreeBase can be
  #  temperamental
  for (i in 1:trys) {
    res <- try (expr, TRUE)
    if (any (class (res) == 'try-error')) {
      error.message <- geterrmessage ()
      if (verbose) {
        cat ('\n.... Connection error occurred:')
        cat (error.message)
        cat (paste0 ('\n.... waiting [', wait,'] seconds'))
      }
      Sys.sleep (wait)
    } else {
      return (res)
    }
  }
  stop ('Unable to connect')
}