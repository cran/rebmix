.legendA <- function(s, text, col, pch, error)
{
  for (i in s:1) {
    legend <- paste(bquote(.(text[1:i])), sep = ""); j <- i

    if (i < s) {
      legend <- c(legend, "..."); j <- i + 1
    }

    if (error) {
      legend <- c(legend, "Error")
    }

    w <- legend("bottom",
      legend = legend,
      col = if (error) c(col[1:j], "black") else col[1:j],
      lty = 0,
      pch = if (error) c(pch[1:j], 1) else pch[1:j],
      bty = "n",
      cex = 1.0,
      y.intersp = 0,
      plot = FALSE,
      horiz = TRUE,
      inset = c(0, 0),
      x.intersp = 0.25,
      xpd = TRUE)$rect$w

    usr <- par("usr")

    if (w <= usr[2] - usr[1]) {
      break
    }
  }

  legend("bottom",
    legend = legend,
    col = if (error) c(col[1:j], "black") else col[1:j],
    lty = 0,
    pch = if (error) c(pch[1:j], 1) else pch[1:j],
    bty = "n",
    cex = 1.0,
    y.intersp = 0,
    horiz = TRUE,
    inset = c(0, 0),
    x.intersp = 0.25,
    xpd = TRUE)
} ## .legendA
