# 16/07/2014
# D.J. Bennett
# Comparing simulated and empirical trees

# LIBS
source (file.path ('tools', 'compare_tools.R'))
source (file.path ('tools', 'misc_tools.R'))

# FUNCTIONS
lmPlot <- function (x_string, y_string, x_lab, y_lab) {
  # Use ggplot to plot relationship with lm
  p <- ggplot (stats, aes_string (x=x_string, y=y_string))
  p <- p + geom_point () + stat_smooth (method='lm') +
    ylab (y_lab) + xlab (x_lab) +
    theme_bw () + theme (text=element_text(size=24))
  print (p)
}

# DIRS
data.dir <- file.path ('data', 'treestats')

# PARAMETERS
analysis.name <- 'analysis_3'  # 1000 trees with different birth-death parameters
res.dir <- file.path ('results', analysis.name)

# INPUT
stats <- readIn (analysis.name)

# PLOT
cor (stats$birth, stats$psv)
lmPlot ('birth', 'psv', 'Birth-Death', 'PSV')
cor (stats$birth, stats$gamma)
lmPlot ('birth', 'gamma', 'Birth-Death', expression (gamma))
cor (stats$eps, stats$gamma)
lmPlot ('eps', 'gamma', expression (epsilon), expression (gamma))
cor (stats$eps, stats$psv)
lmPlot ('eps', 'psv', expression (epsilon), 'PSV')
