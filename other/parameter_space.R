# 16/04/2015
# Parameter space of sig-eps

# QUADRANTS
# empty plot
xlab1 <- expression(epsilon)
xlab2 <- expression(paste('Increasing extinction rate for ED species') %->% '')
ylab1 <- expression(sigma)
ylab2 <- expression(paste('Increasing speciation rate for ED species') %->% '')
tiff("~/Desktop/figure.tiff", width=9, height=9, units="cm",
     res=1200)
par (mar=(c(2, 2, 0, 0) + 1), mgp=c(3, 0.35, 0))
plot (x=0, y=0, pch=19, xlim=c(-1,1),
      ylim=c(-1,1), yaxt="n", xaxt="n", xlab="", ylab="")
# plot labels at different positions
abline (v=0)
abline (h=0)
sig <- c (-0.5, -0.5, 0.5, 0.5)
eps <- c (-0.5, 0.5, -0.5, 0.5)
labels <- c ('Pan', 'PF', 'DE', 'Eph')
text (sig, eps, labels=labels, cex=1)
sig <- c (-1, -1, 1, 1, -1)
eps <- c (-1, 1, -1, 1, 0)
points (eps, sig, pch='X', cex = 0.5)
axis(1, at=c(-1.0, 0, 1.0), cex=0.25, cex.axis=0.5, lwd.ticks=0.5)
axis(2, at=c(-1.0, 0, 1.0), cex=0.25, cex.axis=0.5, lwd.ticks=0.5)
mtext(xlab1, side=1, line=1.2, cex=1)
mtext(xlab2, side=1, line=1.75, cex=0.5)
mtext(ylab1, side=2, line=1.75, cex=1)
mtext(ylab2, side=2, line=1.2, cex=0.5)
dev.off()

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