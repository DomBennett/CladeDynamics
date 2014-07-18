## 18/07/2014
## D.J. Bennett
## Sanity checking code

## Testing the ED strength metric
# Poisson distributed ED values with some noise
eds <- sort (rpois (100, 1) + rnorm (100, mean = 1, sd = 0.1))
# a vector of strengths for calculating probability
ed.strength <- seq (-1, 1, 0.5)
colours <- rainbow (length (ed.strength), alpha = 0.7)
plot (x = c (min (eds), max (eds)), y = c (0, max (eds)/sum(eds)),
      col = NULL, ylab = "Diversification Probabilty",
      xlab = "Evolutionary Distinctiveness")
for (i in 1:length (ed.strength)) {
  probs <- eds^ed.strength[i]
  probs <- probs/sum (probs)
  lines (y = probs, x = eds, col = colours[i], lwd = 2)
}
legend ("top", as.character (ed.strength), col = colours, lwd = 1, cex = 0.5)