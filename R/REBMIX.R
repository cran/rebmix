RNGMIX <- function(Dataset = NULL,
  rseed = -1,
  n = NULL,
  Theta = NULL)
{
  i <- length(Dataset)
    
  InpRNGMIX <- ""

  for (j in 1:i) {
    InpRNGMIX <- paste(InpRNGMIX, "Dataset = ", Dataset[j], "\n", sep = "")
  }

  InpRNGMIX <- paste(InpRNGMIX, "rseed = ", rseed, "\n", sep = "")

  n <- rbind(n)

  nTheta <- rbind(n, Theta)

  d <- as.integer(nrow(Theta) / 3); c <- ncol(Theta)

  nj <- 3 * d + 1

  for (i in 1:c) {
    InpRNGMIX <- paste(InpRNGMIX, "ntheta = [", sep = "")

    for (j in 1:nj) {
      if (j == nj) {
        InpRNGMIX <- paste(InpRNGMIX, nTheta[j, i], "]\n", sep = "")
      }
      else {
        InpRNGMIX <- paste(InpRNGMIX, nTheta[j, i], ", ", sep = "")  
      }
    }
  }

  InpRNGMIX <- paste(InpRNGMIX, "Save = OutRNGMIX.txt\n", sep = "")

  InpRNGMIX <- paste(InpRNGMIX, "Run = RNGMIX\n", sep = "")

  write.table(InpRNGMIX, 
    file = "InpRNGMIX.txt", 
    append = FALSE, 
    quote = FALSE, 
    sep = "\t",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE)

  output <- .C("RRNGMIX",
    file = as.character("InpRNGMIX.txt"), 
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in RNGMIX!"); return(NA)
  }

  RNGMIX <- list(2)

  RNGMIX$w <- data.frame(n / sum(n))
    
  rownames(RNGMIX$w) <- "w"
  colnames(RNGMIX$w) <- paste("comp", 1:c, sep = "")

  RNGMIX$Theta <- data.frame(rbind(Theta))

  if (d == 1) {
    rownames(RNGMIX$Theta) <- c("pdf", "theta1", "theta2")
  } 
  else {
    rownames(RNGMIX$Theta) <- paste(rep(c("pdf", "theta1.", "theta2."), d), rep(1:d, each = 3), sep = "")
  }

  colnames(RNGMIX$Theta) <- paste("comp", 1:c, sep = "")

  rm(list = ls()[!(ls() %in% c("RNGMIX"))])

  class(RNGMIX) <- "RNGMIX"
 
  return(RNGMIX)
} ## RNGMIX

REBMIX <- function(Dataset = NULL, 
  Preprocessing = NULL, 
  D = 0.025, 
  cmax = 15,
  InformationCriterion = "AIC",
  pdf = NULL,
  K = NULL,
  Rmin = 0.001,
  ar = 0.1,
  Restraints = "loose")
{
  i <- length(Dataset)
    
  InpREBMIX <- ""

  for (j in 1:i) {
    InpREBMIX <- paste(InpREBMIX, "Dataset = ", Dataset[j], "\n", sep = "")
  }

  InpREBMIX <- paste(InpREBMIX, "Preprocessing = ", Preprocessing, "\n", sep = "")

  InpREBMIX <- paste(InpREBMIX, "D = ", D, "\n", sep = "")

  InpREBMIX <- paste(InpREBMIX, "cmax = ", cmax, "\n", sep = "")

  InpREBMIX <- paste(InpREBMIX, "InformationCriterion = ", InformationCriterion, "\n", sep = "")

  i <- length(pdf)

  InpREBMIX <- paste(InpREBMIX, "pdf = [", sep = "")

  for (j in 1:i) {
    InpREBMIX <- paste(InpREBMIX, pdf[j], sep = "")

    if (j < i) {
      InpREBMIX <- paste(InpREBMIX, ", ", sep = "")
    } 
    else {
      InpREBMIX <- paste(InpREBMIX, "]\n", sep = "")
    }
  }

  i <- length(K)

  InpREBMIX <- paste(InpREBMIX, "K = [", sep = "") 

  for (j in 1:i) {
    InpREBMIX <- paste(InpREBMIX, K[j], sep = "")

    if (j < i) {
      InpREBMIX <- paste(InpREBMIX, ", ", sep = "")
    } 
    else {
      InpREBMIX <- paste(InpREBMIX, "]\n", sep = "")
    }
  }

  InpREBMIX <- paste(InpREBMIX, "Rmin = ", Rmin, "\n", sep = "")

  InpREBMIX <- paste(InpREBMIX, "ar = ", ar, "\n", sep = "")

  InpREBMIX <- paste(InpREBMIX, "Restraints = ", Restraints, "\n", sep = "")

  InpREBMIX <- paste(InpREBMIX, "Save = OutREBMIX.txt\n", sep = "")

  InpREBMIX <- paste(InpREBMIX, "Run = REBMIX\n", sep = "")

  write.table(InpREBMIX, 
    file = "InpREBMIX.txt", 
    append = FALSE, 
    quote = FALSE, 
    sep = "\t",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE)

  output <- .C("RREBMIX",
    file = as.character("InpREBMIX.txt"), 
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in REBMIX!"); return(NA)
  }

  OutREBMIX_1 <- read.table("OutREBMIX_1.txt", 
    header = TRUE, 
    sep = "\t",
    quote = "",
    dec = ".",
    check.names = FALSE,
    blank.lines.skip = TRUE,
    stringsAsFactors = FALSE)

  OutREBMIX_2 <- read.table("OutREBMIX_2.txt", 
    header = TRUE, 
    sep = "\t",
    quote = "",
    dec = ".",
    check.names = FALSE,
    blank.lines.skip = TRUE,
    stringsAsFactors = FALSE)

  ni <- nrow(OutREBMIX_1); d <- as.integer((ncol(OutREBMIX_2) - 2) / 3)

  REBMIX <- list(3)

  REBMIX$w <- list(ni)
  REBMIX$Theta <- list(ni)
  REBMIX$summary <- OutREBMIX_1

  for (i in 1:ni) {
    TmpREBMIX_2 <- OutREBMIX_2[OutREBMIX_2[ ,1] %in% OutREBMIX_1[i, 1], , drop = FALSE]

    c <- nrow(TmpREBMIX_2)
        
    TmpREBMIX_2 <- TmpREBMIX_2[!(colnames(TmpREBMIX_2) %in% c("Dataset"))]

    TmpREBMIX_2 <- as.data.frame(t(TmpREBMIX_2), stringsAsFactors = FALSE)

    if (c == 1) {
      colnames(TmpREBMIX_2) <- c("comp")
    } 
    else {
      colnames(TmpREBMIX_2) <- paste("comp", 1:c, sep = "")
    }

    REBMIX$w[[i]] <- TmpREBMIX_2[rownames(TmpREBMIX_2) %in% c("w"), , drop = FALSE]
    REBMIX$Theta[[i]] <- TmpREBMIX_2[!(rownames(TmpREBMIX_2) %in% c("w")), , drop = FALSE]
  }

  rm(list = ls()[!(ls() %in% c("REBMIX"))])

  class(REBMIX) <- "REBMIX"
 
  return(REBMIX)
} ## REBMIX 

.densKNearestNeighbour.x <- function(x, k, Rmin, hx) 
{
  output <- .C("RdensKNearestNeighbourX",
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    k = as.integer(k),
    Rmin = as.double(Rmin),
    hx = as.double(hx),
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in densKNearestNeighbour.x!"); return(NA)
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densKNearestNeighbour.x

.densKNearestNeighbour.xy <- function(x, y, k, Rmin, hx, hy) 
{
  output <- .C("RdensKNearestNeighbourXY",
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    k = as.integer(k),
    Rmin = as.double(Rmin),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in densKNearestNeighbour.xy!"); return(NA)
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densKNearestNeighbour.xy

.densParzenWindow.x <- function(x, hx) 
{
  output <- .C("RdensParzenWindowX",
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    hx = as.double(hx),
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in densParzenWindow.x!"); return(NA)
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densParzenWindow.x

.densParzenWindow.xy <- function(x, y, hx, hy) 
{
  output <- .C("RdensParzenWindowXY",
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in densParzenWindow.xy!"); return(NA)
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densParzenWindow.xy

.densHistogram.x <- function(k, x, x0, hx) 
{
  output <- .C("RdensHistogramX",
    k = as.integer(k),
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    x0 = as.double(x0),
    hx = as.double(hx),
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in densHistogram.x!"); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.x 

.densHistogram.xy <- function(k, x, y, x0, y0, hx, hy) 
{
  output <- .C("RdensHistogramXY",
    k = as.integer(k),
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    x0 = as.double(x0),
    y0 = as.double(y0),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    NAOK = FALSE,
    DUP = TRUE)

  if (output$error == 1) {
    message("Error in densHistogram.xy!"); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k
  length(output$z) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.xy

.dfmix.x <- function(x, w, xpdf, xtheta0, xtheta1, ...) 
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    C <- toupper(xpdf[i])

    if (C == "NORMAL") {
      fix <- dnorm(x, mean = xtheta0[i], sd = xtheta1[i], ...)
    } 
    else 
    if (C == "LOGNORMAL") {
      fix <- dlnorm(x, meanlog = xtheta0[i], sdlog = xtheta1[i], ...)
    } 
    else 
    if (C == "WEIBULL") {
      fix <- dweibull(x, scale = xtheta0[i], shape = xtheta1[i], ...)
    } 
    else 
    if (C == "BINOMIAL") {
      fix <- dbinom(x, size = as.integer(xtheta0[i]), prob = xtheta1[i], ...)
    }

    f <- f + w[i] * fix
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## .dfmix.x

.dfmix.xy <- function(x, y, w, xpdf, xtheta0, xtheta1, ypdf, ytheta0, ytheta1, ...) 
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    C <- toupper(xpdf[i])

    if (C == "NORMAL") {
      fix <- dnorm(x, mean = xtheta0[i], sd = xtheta1[i], ...)
    } 
    else 
    if (C == "LOGNORMAL") {
      fix <- dlnorm(x, meanlog = xtheta0[i], sdlog = xtheta1[i], ...)
    } 
    else 
    if (C == "WEIBULL") {
      fix <- dweibull(x, scale = xtheta0[i], shape = xtheta1[i], ...)
    } 
    else 
    if (C == "BINOMIAL") {
      fix <- dbinom(x, size = as.integer(xtheta0[i]), prob = xtheta1[i], ...)
    }

    C <- toupper(ypdf[i])

    if (C == "NORMAL") {
      fiy <- dnorm(y, mean = ytheta0[i], sd = ytheta1[i], ...)
    } 
    else 
    if (C == "LOGNORMAL") {
      fiy <- dlnorm(y, meanlog = ytheta0[i], sdlog = ytheta1[i], ...)
    } 
    else 
    if (C == "WEIBULL") {
      fiy <- dweibull(y, scale = ytheta0[i], shape = ytheta1[i], ...)
    } 
    else 
    if (C == "BINOMIAL") {
      fiy <- dbinom(y, size = as.integer(ytheta0[i]), prob = ytheta1[i], ...)
    }

    f <- f + w[i] * fix * fiy
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## .dfmix.xy

plot.REBMIX <- function(x,
  pos = 1,
  nrow = 2, 
  ncol = 2, 
  npts = 200,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19,
  plot.type = "p",
  contour.drawlabels = FALSE,
  contour.labcex = 0.8, 
  contour.method = "flattest",
  contour.nlevels = 12, ...) 
{
  ni <- ncol(x$summary[ , , drop = FALSE])

  d <- as.integer(nrow(x$Theta[[pos]]) / 3)

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)
    
  n <- d * (d - 1) / 2

  if (n == 0) {
    ncol = 1; nrow = 1
  }

  opar <- par(mfrow = c(nrow, ncol),
    cex = cex, 
    cex.axis = 1.0,
    fg = fg, 
    lty = lty,
    lwd = lwd,
    mar = c(1.2, 1.2, 1.2, 1.2),
    oma = c(1.2, 0.2, 0.2, 0.2),
    pty = pty, 
    tcl = tcl, ...)

  if (.Device == "tikz output") {
    item <- list(1)

    item[[1]] <- "$\\textrm{Dataset}$"
    item[[2]] <- "$\\; = \\;$"
    item[[3]] <- paste("$\\textrm{", x$summary[pos, 1, drop = FALSE], "}$, ", sep = "")

    item[[4]] <- "$\\textrm{Preprocessing}$"
    item[[5]] <- "$\\; = \\;$"
    item[[6]] <- paste("$\\textrm{", x$summary[pos, 2, drop = FALSE], "}$, ", sep = "")

    item[[7]] <- "$\\textrm{Restraints}$"
    item[[8]] <- "$\\; = \\;$"
    item[[9]] <- paste("$\\textrm{", x$summary[pos, 7, drop = FALSE], "}$, ", sep = "")

    item[[10]] <- "$D$"
    item[[11]] <- "$\\; = \\;$"
    item[[12]] <- paste("$", x$summary[pos, 3, drop = FALSE], "$, ", sep = "")

    item[[13]] <- "$c_{\\mathrm{max}}$"
    item[[14]] <- "$\\; = \\;$"
    item[[15]] <- paste("$", x$summary[pos, 4, drop = FALSE], "$, ", sep = "")

    item[[16]] <- "$a_{\\mathrm{r}}$"
    item[[17]] <- "$\\; = \\;$"
    item[[18]] <- paste("$", x$summary[pos, 6, drop = FALSE], "$, ", sep = "")

    item[[19]] <- "$c$"
    item[[20]] <- "$\\; = \\;$"
    item[[21]] <- paste("$", x$summary[pos, 8, drop = FALSE], "$, ", sep = "")

    item[[22]] <- "$k$"
    item[[23]] <- "$\\; = \\;$"
    item[[24]] <- paste("$", x$summary[pos, 9, drop = FALSE], "$, ", sep = "")

    item[[25]] <- "$t_{\\mathrm{c}}$"
    item[[26]] <- "$\\; = \\;$"
    item[[27]] <- paste("$", x$summary[pos, ni - 2, drop = FALSE], "$, ", sep = "")

    item[[28]] <- paste("$\\mathrm{", x$summary[pos, 5, drop = FALSE], "}$", sep = "")
    item[[29]] <- "$\\; = \\;$"
    item[[30]] <- paste("$", format(x$summary[pos, ni - 1, drop = FALSE], digits = 3, trunc = TRUE), "$, ", sep = "")

    item[[31]] <- "$\\mathrm{log}\\, L$"
    item[[32]] <- "$\\; = \\;$"
    item[[33]] <- paste("$", format(x$summary[pos, ni, drop = FALSE], digits = 3, trunc = TRUE), "$.", sep = "")

    i <- 1; legend <- list(1); legend[[i]] <- item[[1]]

    for (j in 2:33) {
      legendwidth <- strwidth(paste(legend[[i]], item[[j]], sep = ""), units = "figure", cex = 1.0)
  
      if (legendwidth > ncol) {
        i <- i + 1; legend[[i]] <- item[[j]]
      }
      else {
        legend[[i]] <- paste(legend[[i]], item[[j]], sep = "")
      }
    }
  }
  else {
    item <- list(1)

    item[[1]] <- "Dataset"
    item[[2]] <- " = "
    item[[3]] <- paste(x$summary[pos, 1, drop = FALSE], ", ", sep = "")

    item[[4]] <- "Preprocessing"
    item[[5]] <- " = "
    item[[6]] <- paste(x$summary[pos, 2, drop = FALSE], ", ", sep = "")

    item[[7]] <- "Restraints"
    item[[8]] <- " = "
    item[[9]] <- paste(x$summary[pos, 7, drop = FALSE], ", ", sep = "")

    item[[10]] <- "D"
    item[[11]] <- " = "
    item[[12]] <- paste(x$summary[pos, 3, drop = FALSE], ", ", sep = "")

    item[[13]] <- bquote(c[max])
    item[[14]] <- " = "
    item[[15]] <- paste(x$summary[pos, 4, drop = FALSE], ", ", sep = "")

    item[[16]] <- bquote(a[r])
    item[[17]] <- " = "
    item[[18]] <- paste(x$summary[pos, 6, drop = FALSE], ", ", sep = "")

    item[[19]] <- "c"
    item[[20]] <- " = "
    item[[21]] <- paste(x$summary[pos, 8, drop = FALSE], ", ", sep = "")

    item[[22]] <- "k"
    item[[23]] <- " = "
    item[[24]] <- paste(x$summary[pos, 9, drop = FALSE], ", ", sep = "")

    item[[25]] <- bquote(t[c])
    item[[26]] <- " = "
    item[[27]] <- paste(x$summary[pos, ni - 2, drop = FALSE], ", ", sep = "")

    item[[28]] <- as.character(x$summary[pos, 5, drop = FALSE])
    item[[29]] <- " = "
    item[[30]] <- paste(format(x$summary[pos, ni - 1, drop = FALSE], digits = 3, trunc = TRUE), ", ", sep = "")

    item[[31]] <- "log L"
    item[[32]] <- " = "
    item[[33]] <- paste(format(x$summary[pos, ni, drop = FALSE], digits = 3, trunc = TRUE), ".", sep = "")

    i <- 1; legend <- list(1); legend[[i]] <- bquote(.(item[[1]]))

    for (j in 2:33) {
      legendwidth <- strwidth(bquote(paste(.(legend[[i]]), .(item[[j]]), sep = "")), units = "figure", cex = 1.0)
  
      if (legendwidth > ncol) {
        i <- i + 1; legend[[i]] <- item[[j]]
      }
      else {
        legend[[i]] <- bquote(paste(.(legend[[i]]), .(item[[j]]), sep = ""))
      }
    }
  }

  par(oma = c(length(legend) + 0.2, 0.2, 0.2, 0.2))

  ey <- read.table(paste(x$summary[pos, 1, drop = FALSE], ".txt", sep = ""), 
    header = FALSE, 
    sep = "\t",
    quote = "",
    dec = ".",
    check.names = FALSE,
    blank.lines.skip = TRUE)

  y0 <- array(data = 0.0, dim = d, dimnames = NULL)
  h <- array(data = 0.0, dim = d, dimnames = NULL)

  lim <- array(data = 0.0, dim = c(2, d), dimnames = NULL)

  py <- array(data = 0.0, dim = c(npts / ncol, d), dimnames = NULL)

  C <- toupper(x$summary[pos, 2, drop = FALSE])

  for (i in 1:d) {
    if (C == "HISTOGRAM") {
      k <- as.numeric(x$summary[pos, 9, drop = FALSE])
      y0[i] <- as.numeric(x$summary[pos, 9 + i, drop = FALSE])
      h[i] <- as.numeric(x$summary[pos, 9 + d + i, drop = FALSE])

      lim[, i] <- range(ey[, i], finite = TRUE) + 0.5 * h[i]
    } 
    else 
    if (C == "PARZEN WINDOW") {
      h[i] <- as.numeric(x$summary[pos, 9 + i, drop = FALSE])

      lim[, i] <- range(ey[, i], finite = TRUE)
    } 
    else
    if (C == "K-NEAREST NEIGHBOUR") {
      k <- as.numeric(x$summary[pos, 9, drop = FALSE])
      Rmin <- as.numeric(x$summary[pos, 10, drop = FALSE])

      h[i] <- as.numeric(x$summary[pos, 10 + i, drop = FALSE])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }

    py[, i] <- seq(from = lim[1, i], to = lim[2, i], length.out = length(py[, i])) 
  }

  x$w[[pos]] <- as.numeric(x$w[[pos]])

  if (n > 0) {
    ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
      space = "rgb",
      interpolate = "linear")

    figno <- 0

    for (i in 1:(d - 1)) {
      xpdf <- as.character(x$Theta[[pos]][3 * i - 2, ])
      xtheta0 <- as.numeric(x$Theta[[pos]][3 * i - 1, ])
      xtheta1 <- as.numeric(x$Theta[[pos]][3 * i, ])

      for (j in (i + 1):d) {
        if (C == "HISTOGRAM") {
          edens <- .densHistogram.xy(k, ey[, i], ey[, j], y0[i], y0[j], h[i], h[j]) 
        } 
        else 
        if (C == "PARZEN WINDOW") {
          edens <- .densParzenWindow.xy(ey[, i], ey[, j], h[i], h[j])
        } 
        else
        if (C == "K-NEAREST NEIGHBOUR") {
          edens <- .densKNearestNeighbour.xy(ey[, i], ey[, j], k, Rmin, h[i], h[j])
        }

        ypdf <- as.character(x$Theta[[pos]][3 * j - 2, ])
        ytheta0 <- as.numeric(x$Theta[[pos]][3 * j - 1, ])
        ytheta1 <- as.numeric(x$Theta[[pos]][3 * j, ])

        pdens <- outer(py[, i], py[, j], ".dfmix.xy", x$w[[pos]], xpdf, xtheta0, xtheta1, ypdf, ytheta0, ytheta1)

        zlim <- range(edens$z, finite = TRUE); zmax <- max(zlim[2], pdens)

        zlim <- zlim / zmax

        plot(x = edens$x, 
          y = edens$y, 
          type = plot.type,
          main = "", 
          sub = "",   
          xlab = "",
          ylab = "",  
          col = rgb(ramp(edens$z / zmax), maxColorValue = 255),
          axes = FALSE,
          cex = plot.cex,
          pch = plot.pch)

        levels <- 10^seq(from = log(zlim[1]), to = log(zlim[2]), length.out = contour.nlevels)

        contour(x = py[, i],
          y = py[, j],
          z = pdens / zmax, 
          levels = levels,
          xlim = lim[, i],
          ylim = lim[, j],
          zlim = zlim,
          labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
          axes = FALSE, frame.plot = FALSE,
          col = rgb(ramp(levels), maxColorValue = 255), 
          add = TRUE)

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

        if (.Device == "tikz output") {
          text <- paste("$y_{", i, "} - y_{", j, "}$", sep = "")
        }
        else {
          text <- bquote(y[.(i)] - y[.(j)])
        }

        mtext(text = text, 
          side = 1, 
          line = 0, 
          outer = FALSE,
          adj = 0.5,
          padj = 0.0,
          cex = cex)

        figno <- figno + 1

        if ((figno == nrow * ncol) || ((i == d - 1) && (j == d))) {
          for (k in 1:length(legend)) {
            mtext(text = legend[[k]],
              side = 1, 
              line = k - 1, 
              outer = TRUE,
              adj = 0.5,
              padj = 0.0,
              cex = cex)
          }

          figno <- 0
        }
      }
    }
  } 
  else {
    if (C == "HISTOGRAM") {
      edens <- .densHistogram.x(k, ey[, 1], y0[1], h[1])
    } 
    else 
    if (C == "PARZEN WINDOW") {
      edens <- .densParzenWindow.x(ey[, 1], h[1])
    } 
    else
    if (C == "K-NEAREST NEIGHBOUR") {
      edens <- .densKNearestNeighbour.x(ey[, 1], k, Rmin, h[1])
    }

    xpdf <- as.character(x$Theta[[pos]][1, ])
    xtheta0 <- as.numeric(x$Theta[[pos]][2, ])
    xtheta1 <- as.numeric(x$Theta[[pos]][3, ])

    pdens <- .dfmix.x(py[, 1], x$w[[pos]], xpdf, xtheta0, xtheta1)
        
    ylim <- c(0.0, max(edens$y, pdens))

    plot(x = edens$x, 
      y = edens$y, 
      type = plot.type,
      main = "", 
      sub = "",   
      xlab = "",
      ylab = "",  
      ylim = ylim,
      col = "black",
      axes = FALSE,
      cex = plot.cex,
      pch = plot.pch)

    points(x = py[, 1], 
      y = pdens, 
      type = "l", 
      col = "black")

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

    if (.Device == "tikz output") {
      text <- "$y_{1} - f(y_{1})$" 
    }
    else {
      text <- bquote(y[1] - f(y[1]))
    }

    mtext(text = text, 
      side = 1, 
      line = 0, 
      outer = FALSE,
      adj = 0.5,
      padj = 0.0,
      cex = cex)

    for (i in 1:length(legend)) {
      mtext(text = legend[[i]],
        side = 1, 
        line = i - 1, 
        outer = TRUE,
        adj = 0.5,
        padj = 0.0,
        cex = cex)
    }
  }

  par(opar)

  rm(list = ls())
} ## plot.REBMIX
