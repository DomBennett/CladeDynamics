## 12/09/2014
## D.J. Bennett
## Misc tools

closeDevices <- function () {
  # make sure all graphical devices are off
  while (!is.null (dev.list ())) {
    dev.off ()
  }
}