# 16/04/2015
# Parameter space of sig-eps

# empty plot
par (mar=(c(4, 4, 4, 4) + 1))
xlab <- expression (paste (epsilon, ' (Increasing extinction rate for ED species -->)'))
ylab <- expression (paste (sigma, ' (Increasing speciation rate for ED species -->)'))
plot (x=0, y=0, pch=19, xlab=xlab, ylab=ylab, xlim=c(-1,1),
      ylim=c(-1,1), cex.lab=1.5)
# plot labels at different positions
abline (v=0)
abline (h=0)
sig <- c (-0.5, -0.5, 0.5, 0.5)
eps <- c (-0.5, 0.5, -0.5, 0.5)
labels <- c ('Pan', 'PF', 'DE', 'Eph')
text (sig, eps, labels=labels, cex=2)
