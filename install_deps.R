# 30/06/2015
# Dom Bennett
# Ensure deps are installed

# START
user.input <- readline('WARNING: this script will automaticaly install packages.
                       Are you sure you want to continue? (y/n) ')
if (user.input != 'y') {
  stop ('Execution halted')
}

# GET INSTALLED PACKAGES
packages <- installed.packages ()[ ,1]

# CHECK AND INSTALL
counter <- 0
if (!'devtools' %in% packages) {
  install.packages ('devtools')
  counter <- counter + 1
}
if (!'MoreTreeTools' %in% packages) {
  library (devtools)
  install_github ('https://github.com/DomBennett/MoreTreeTools.git')
  counter <- counter + 1
}
if (!'ape' %in% packages) {
  install.packages('ape')
  counter <- counter + 1
}
if (!'plyr' %in% packages) {
  install.packages('plyr')
  counter <- counter + 1
}
if (!'geiger' %in% packages) {
  install.packages('geiger')
  counter <- counter + 1
}
if (!'apTreeshape' %in% packages) {
  install.packages('apTreeshape')
  counter <- counter + 1
}
if (!'caper' %in% packages) {
  install.packages('caper')
  counter <- counter + 1
}
if (!'ggplot2' %in% packages) {
  install.packages('ggplot2')
  counter <- counter + 1
}
if (!'foreach' %in% packages) {
  install.packages('foreach')
  counter <- counter + 1
}
if (!'doMC' %in% packages) {
  install.packages('doMC')
  counter <- counter + 1
}
if (!'outliers' %in% packages) {
  install.packages('outliers')
  counter <- counter + 1
}
if (!'treebase' %in% packages) {
  install.packages('treebase')
  counter <- counter + 1
}

# END
cat ('\nComplete! Installed [', counter, '] packages.', sep = '')