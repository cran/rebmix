setMethod("plot",
          signature(x = "RCLRMIX", y = "missing"),
function(x,
  y,
  s = expression(c),
  nrow = 1,
  ncol = 1,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLRMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(nrow)) {
    stop(sQuote("nrow"), " integer is requested!", call. = FALSE)
  }

  if (nrow < 1) {
    stop(sQuote("nrow"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(ncol)) {
    stop(sQuote("ncol"), " integer is requested!", call. = FALSE)
  }

  if (ncol < 1) {
    stop(sQuote("ncol"), " must be greater than 0!", call. = FALSE)
  }

  d <- length(x@x@Variables)

  Zp <- as.numeric(levels(x@Zp))[x@Zp]
  Zt <- as.numeric(levels(x@Zt))[x@Zt]

  zlim <- c(0, max(1, max(Zp) - 1)); zmax <- zlim[2]

  c <- x@c; s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }
  
  unique.Zp <- unique(Zp)
  
  set <- which(x@from %in% unique.Zp)

  from <- x@from[set]; to <- x@to[set]
  
  i <- length(unique.Zp)

  while (i > s) {
    i <- i - 1
    
    Zp[Zp == from[i]] <- to[i]
  }

  sort.unique.Zp <- sort(unique.Zp)

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2

  opar <- list(); ipar <- 1
  
  opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1  

  par(mfrow = c(nrow, ncol),
    cex = cex,
    cex.axis = 1.0,
    fg = fg,
    lty = lty,
    lwd = lwd,
    mar = c(1.2, 1.2, 1.2, 1.2),
    oma = c(1.2, 0.2, 0.2, 0.2),
    pty = pty,
    tcl = tcl, ...)

  par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

  if (as.character(class(x@Dataset)) == "data.frame") {
    ey <- as.matrix(x@Dataset)
  }
  else
  if (as.character(class(x@Dataset)) == "Histogram") {
    ey <- as.matrix(x@Dataset@Y[, 1:d])
  }  
  
  ep <- Zp - 1

  error <- is.error(Zt, Zp)

  ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
    space = "rgb",
    interpolate = "linear")

  plot.col <- rgb(ramp(ep / zmax), maxColorValue = 255)

  if (i > 1) {
    plot.mul <- ifelse((Zp == from[i - 1]) | (Zp == to[i - 1]), 1.5, 1.0)
  }
  else {
    plot.mul <- rep(1.0, length(Zp))
  }

  legend.col <- rgb(ramp((sort.unique.Zp - 1) / zmax), maxColorValue = 255)

  legend.text <- as.character(sort.unique.Zp)

  legend.pch <- rep(plot.pch, s)

  which.error <- which(error == 1); which.not.error <- which(error != 1)

  if (N > 0) {
    figno <- 0

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[which.not.error, i],
          y = ey[which.not.error, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          xlim = range(ey[, i]),
          ylim = range(ey[, j]),
          col = plot.col[which.not.error],
          axes = FALSE,
          lwd = 1,
          cex = plot.cex * plot.mul[which.not.error],
          pch = plot.pch)

        points(x = ey[which.error, i],
          y = ey[which.error, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          col = "black",
          lwd = 0.6,
          cex = 0.6 * plot.cex * plot.mul[which.error],
          pch = 1)

        box(col = fg, lty = "solid", lwd = 1)

        axis(side = 3,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        axis(side = 2,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        text <- bquote(y[.(i)] - y[.(j)])

        mtext(text = text,
          side = 1,
          line = 0,
          outer = FALSE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)

        figno <- figno + 1

        if ((figno == nrow * ncol) || ((i == d - 1) && (j == d))) {
          par(fig = c(0, 1, 0, 1),
            oma = c(0, 0, 0, 0),
            mar = c(0, 0, 0, 0),
            new = TRUE)

          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

          .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = sum(error) != 0)

          par(mfrow = c(nrow, ncol),
            cex = cex,
            cex.axis = 1.0,
            fg = fg,
            lty = lty,
            lwd = lwd,
            mar = c(1.2, 1.2, 1.2, 1.2),
            oma = c(1.2, 0.2, 0.2, 0.2),
            pty = pty,
            tcl = tcl, ...)

          par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }

        opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
      }
    }
  }
  else {
    plot(x = ey[which.not.error, 1],
      y = Zp[which.not.error],
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      xlim = range(ey[, 1]),
      ylim = range(Zp),
      col = plot.col[which.not.error],
      axes = FALSE,
      lwd = 1,
      cex = plot.cex * plot.mul[which.not.error],
      pch = plot.pch)

    points(x = ey[which.error, 1],
      y = Zp[which.error],
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      col = "black",
      lwd = 0.6,
      cex = 0.6 * plot.cex * plot.mul[which.error],
      pch = 1)

    box(col = fg, lty = "solid", lwd = 1)

    axis(side = 3,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    axis(side = 2,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    text <- bquote(y[1] - Z[p](y[1]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    par(fig = c(0, 1, 0, 1),
      oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0),
      new = TRUE)

    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

   .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = sum(error) != 0)

    par(mfrow = c(nrow, ncol),
      cex = cex,
      cex.axis = 1.0,
      fg = fg,
      lty = lty,
      lwd = lwd,
      mar = c(1.2, 1.2, 1.2, 1.2),
      oma = c(1.2, 0.2, 0.2, 0.2),
      pty = pty,
      tcl = tcl, ...)

    par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot

setMethod("plot",
          signature(x = "RCLRMVNORM", y = "missing"),
function(x,
  y,
  s = expression(c),
  nrow = 1,
  ncol = 1,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLRMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(nrow)) {
    stop(sQuote("nrow"), " integer is requested!", call. = FALSE)
  }

  if (nrow < 1) {
    stop(sQuote("nrow"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(ncol)) {
    stop(sQuote("ncol"), " integer is requested!", call. = FALSE)
  }

  if (ncol < 1) {
    stop(sQuote("ncol"), " must be greater than 0!", call. = FALSE)
  }

  d <- length(x@x@Variables)

  Zp <- as.numeric(levels(x@Zp))[x@Zp]
  Zt <- as.numeric(levels(x@Zt))[x@Zt]

  zlim <- c(0, max(1, max(Zp) - 1)); zmax <- zlim[2]

  c <- x@c; s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }

  unique.Zp <- unique(Zp)
  
  set <- which(x@from %in% unique.Zp)

  from <- x@from[set]; to <- x@to[set]
  
  i <- length(unique.Zp)

  while (i > s) {
    i <- i - 1
    
    Zp[Zp == from[i]] <- to[i]
  }

  sort.unique.Zp <- sort(unique.Zp)

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2

  opar <- list(); ipar <- 1
  
  opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1

  par(mfrow = c(nrow, ncol),
    cex = cex,
    cex.axis = 1.0,
    fg = fg,
    lty = lty,
    lwd = lwd,
    mar = c(1.2, 1.2, 1.2, 1.2),
    oma = c(1.2, 0.2, 0.2, 0.2),
    pty = pty,
    tcl = tcl, ...)

  par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

  if (as.character(class(x@Dataset)) == "data.frame") {
    ey <- as.matrix(x@Dataset)
  }
  else
  if (as.character(class(x@Dataset)) == "Histogram") {
    ey <- as.matrix(x@Dataset@Y[, 1:d])
  }  

  ep <- Zp - 1

  error <- is.error(Zt, Zp)

  ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
    space = "rgb",
    interpolate = "linear")

  plot.col <- rgb(ramp(ep / zmax), maxColorValue = 255)

  if (i > 1) {
    plot.mul <- ifelse((Zp == from[i - 1]) | (Zp == to[i - 1]), 1.5, 1.0)
  }
  else {
    plot.mul <- rep(1.0, length(Zp))
  }

  legend.col <- rgb(ramp((sort.unique.Zp - 1) / zmax), maxColorValue = 255)

  legend.text <- as.character(sort.unique.Zp)

  legend.pch <- rep(plot.pch, s)

  which.error <- which(error == 1); which.not.error <- which(error != 1)

  if (N > 0) {
    figno <- 0

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[which.not.error, i],
          y = ey[which.not.error, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          xlim = range(ey[, i]),
          ylim = range(ey[, j]),
          col = plot.col[which.not.error],
          axes = FALSE,
          lwd = 1,
          cex = plot.cex * plot.mul[which.not.error],
          pch = plot.pch)

        points(x = ey[which.error, i],
          y = ey[which.error, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          col = "black",
          lwd = 0.6,
          cex = 0.6 * plot.cex * plot.mul[which.error],
          pch = 1)

        box(col = fg, lty = "solid", lwd = 1)

        axis(side = 3,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        axis(side = 2,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        text <- bquote(y[.(i)] - y[.(j)])

        mtext(text = text,
          side = 1,
          line = 0,
          outer = FALSE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)

        figno <- figno + 1

        if ((figno == nrow * ncol) || ((i == d - 1) && (j == d))) {
          par(fig = c(0, 1, 0, 1),
            oma = c(0, 0, 0, 0),
            mar = c(0, 0, 0, 0),
            new = TRUE)

          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

          .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = sum(error) != 0)

          par(mfrow = c(nrow, ncol),
            cex = cex,
            cex.axis = 1.0,
            fg = fg,
            lty = lty,
            lwd = lwd,
            mar = c(1.2, 1.2, 1.2, 1.2),
            oma = c(1.2, 0.2, 0.2, 0.2),
            pty = pty,
            tcl = tcl, ...)

          par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }

        opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
      }
    }
  }
  else {
    plot(x = ey[which.not.error, 1],
      y = Zp[which.not.error],
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      xlim = range(ey[, 1]),
      ylim = range(Zp),
      col = plot.col[which.not.error],
      axes = FALSE,
      lwd = 1,
      cex = plot.cex * plot.mul[which.not.error],
      pch = plot.pch)

    points(x = ey[which.error, 1],
      y = Zp[which.error],
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      col = "black",
      lwd = 0.6,
      cex = 0.6 * plot.cex * plot.mul[which.error],
      pch = 1)

    box(col = fg, lty = "solid", lwd = 1)

    axis(side = 3,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    axis(side = 2,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    text <- bquote(y[1] - Z[p](y[1]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    par(fig = c(0, 1, 0, 1),
      oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0),
      new = TRUE)

    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

   .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = sum(error) != 0)

    par(mfrow = c(nrow, ncol),
      cex = cex,
      cex.axis = 1.0,
      fg = fg,
      lty = lty,
      lwd = lwd,
      mar = c(1.2, 1.2, 1.2, 1.2),
      oma = c(1.2, 0.2, 0.2, 0.2),
      pty = pty,
      tcl = tcl, ...)

    par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot
