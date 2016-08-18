library(ggplot2)

ed_vals <- c(rep((1/6 + 1/4 + 1/2 + 1), 4),
             rep((1/6 + 1/2 + 2), 2), 4)

sigs <- c(-1, -0.5, 0, 0.5, 1)
res <- NULL
for(i in 1:length(sigs)) {
  mdfd <- ed_vals^sigs[i]
  res <- c(res, mdfd/(sum(mdfd)))
}
res <- data.frame(p=res, ed=rep(ed_vals, length(sigs)),
                  sig=rep(sigs, each=length(ed_vals)))
p <- ggplot(res, aes(x=ed, y=p, group=sig, colour=sig)) +
  geom_line(size= 2) + theme_bw() +
  scale_colour_gradient(expression(sigma), low='#ff9999', high="#e60000") +
  xlab("Tip Evolutionary Distinctness") + ylab("P(S)") +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=6),
        legend.text=element_text(size=6), legend.position="top")
tiff("~/Desktop/figure.tiff", width=7.5, height=8, units="cm",
     res=1200)
print(p)
dev.off()


