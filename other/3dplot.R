library('scatterplot3d')

psi <- c (1, -1, -1, 1, 1, -1, -1, 1, 0)
sig <- c (1, -1, 1, -1, 1, -1, 1, -1, 0)
epi <- c (-1, 1, -1, 1, 1, -1, 1, -1, 0)
labels <- c ('PF', 'PF', 'DE', 'DE',
             'Eph', 'Eph', 'Pan', 'Pan',
             'ER')
s3d <- scatterplot3d (x=psi, y=epi, z=sig, pch=19, xlab=expression (psi),
                      ylab=expression (epsilon), zlab=expression (sigma))
s3d.coords <- s3d$xyz.convert(psi, sig, epi)
text(s3d.coords$x, s3d.coords$y, labels=labels,
     cex=.75, pos=4) 
