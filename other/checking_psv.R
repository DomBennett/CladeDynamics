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
cor.test (stats$birth, stats$psv)
lmPlot ('birth', 'psv', 'Birth-Death', 'PSV')
cor.test (stats$birth, stats$gamma)
lmPlot ('birth', 'gamma', 'Birth-Death', expression (gamma))
cor.test (stats$eps, stats$gamma)
lmPlot ('eps', 'gamma', expression (epsilon), expression (gamma))
cor.test (stats$eps, stats$psv)
lmPlot ('eps', 'psv', expression (epsilon), 'PSV')

candidates <- c ("psv", "age", "pd", "ntaxa", "colless", "sackin", 'birth')
for (each in candidates) {
  m <- lm (stats$gamma~stats[ ,each])
  res <- cor (stats$eps, m$residuals)
  cat ('\n', each, ':', res)
}

candidates <- c ("psv", "age", "pd", "ntaxa", "colless", "sackin", 'birth', 'gamma')
for (each in candidates) {
  m <- lm (stats$birth~stats[ ,each])
  res <- summary (m)['coefficients'][[1]][2,'Estimate']
  cat ('\n', each, ':', res)
}

m1 <- lm (stats$colless~stats$gamma)
cor (stats$birth, m1$residuals)
m2 <- lm (stats$psv~m1$residuals)
cor (stats$eps, m1$residuals)
cor (stats$eps, m2$residuals)
cor (stats$eps, stats$psv)
m3 <- lm (m2$residuals ~stats$eps)
summary(m3)
plot (m1$residuals ~stats$eps)

m4 <- lm (stats$psv ~ stats$eps)
summary (m4)
m5 <- lm (stats$gamma ~ stats$eps)
summary (m5)

m6 <- lm (stats$gamma ~ stats$birth)
summary (m6)
m7 <- lm (stats$eps ~ m6$residuals)
summary (m7)
plot (m6$residuals ~ stats$eps)

m8 <- lm (stats$colless~stats$gamma)
plot (stats$birth ~ m1$fitted.values)
m9 <- lm (stats$eps~m1$residuals)
