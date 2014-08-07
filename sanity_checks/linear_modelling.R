## 06/08/2014
## D.J. Bennett
## Working on getting a linear model and find corresponding

## Functions
drawCorresPoints <- function (x, y, predicted, actual, ...) {
  ## Draw lines on a graph showing the corresponding values of
  ##  y given values of x
  # find the minimum values of x and y in plotted space
  min.y <- floor (min (y)) + floor (min (y))
  min.x <- floor (min (x)) + floor (min (x))
  # draw lines
  lines (x = c (actual, actual, min.y), y = c (min.x, predicted, predicted), ...)
}

# Genrate some fake but probable data
balance <- rep (seq (2, -1, -0.1), each = 10)
balance <- balance + rnorm (length (balance), 0, 0.5)
strength <- rep (seq (-1, 2, 0.1), each = 10)
strength <- strength + rnorm (length (strength), 0, 0.5)
plot (strength ~ balance, pch = 19, col = rainbow (3, alpha = 0.7)[3])

# Add linear model line
model <- lm (balance ~ strength)
abline (model)

# Draw on predicted strength values given known nat balance stats
natural.balances <- rnorm (100, 1.5, 0.1)
# draw lower 5%
drawCorresPoints (strength, balance, value = quantile (natural.balances, 0.05),
                  col = 'red', lty = 2)
drawCorresPoints (strength, balance, value = quantile (natural.balances, 0.95),
                  col = 'red', lty = 2)
drawCorresPoints (strength, balance, value = mean (natural.balances),
                  lwd = 2)

drawCorresPoints <- function (model, distribution) {
  ## Plot mean, mean and 5% and 95% quantiles of an x distribution's
  ##  equivalent y values
  .draw <- function (actual, predicted, ...) {
    lines (x = c (actual, actual, min.x),
           y = c (min.y, predicted, predicted), ...)
  }
  # find the minimum values of x and y in plotted space
  min.y <- floor (min (model$model[ ,'y'])) + floor (min (model$model[ ,'y']))
  min.x <- floor (min (model$model[ ,'x'])) + floor (min (model$model[ ,'x']))
  # get actual
  q1.x <- quantile (distribution, 0.05)
  q2.x <- quantile (distribution, 0.95)
  m.x <- mean (distribution)
  # get predicted values
  q1.y <- predict (model, data.frame (x = q1.x))
  q2.y <- predict (model, data.frame (x = q2.x))
  m.y <- predict (model, data.frame (x = m.x))
  # draw lines
  .draw (q1.x, q1.y, col = 'red', lty = 2)
  .draw (q2.x, q2.y, col = 'red', lty = 2)
  .draw (m.x, m.y, lwd = 2)
}