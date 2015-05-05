# 16/04/2015
# Parameter space of sig-eps

# QUADRANTS
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

# DECSEXTANTS
# empty plot
par (mar=(c(4, 4, 4, 4) + 1))
xlab <- expression (paste (epsilon, ' (Increasing extinction rate for ED species') %->% ')')
ylab <- expression (paste (sigma, ' (Increasing speciation rate for ED species') %->% ')')
plot (x=0, y=0, pch=19, xlab=xlab, ylab=ylab, xlim=c(-1,1),
      ylim=c(-1,1), cex.lab=1.5)
# plot labels at different positions
abline (v=0, lwd=2)
abline (v=0.5)
abline (v=-0.5)
abline (h=0, lwd=2)
abline (h=0.5)
abline (h=-0.5)
# Pan section
sig <- c (-0.75, -0.75, -0.25, -0.25)
eps <- c (-0.75, -0.25, -0.75, -0.25)
labels <- c ('00', '01', '10', '11')
text (sig, eps, labels=labels, cex=1)
text (-0.5, -0.5, labels='Pan', cex=2)
# PF
sig <- c (-0.75, -0.75, -0.25, -0.25)
eps <- c (0.75, 0.25, 0.75, 0.25)
labels <- c ('01', '00', '11', '10')
text (sig, eps, labels=labels, cex=1)
text (-0.5, 0.5, labels='PF', cex=2)
# DE
sig <- c (0.75, 0.75, 0.25, 0.25)
eps <- c (-0.75, -0.25, -0.75, -0.25)
labels <- c ('10', '11', '00', '01')
text (sig, eps, labels=labels, cex=1)
text (0.5, -0.5, labels='DE', cex=2)
# Eph
sig <- c (0.75, 0.75, 0.25, 0.25)
eps <- c (0.75, 0.25, 0.75, 0.25)
labels <- c ('11', '10', '01', '00')
text (sig, eps, labels=labels, cex=1)
text (0.5, 0.5, labels='Eph', cex=2, bg='white')