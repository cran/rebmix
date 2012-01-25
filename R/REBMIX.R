is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) 
{ 
  warn <- getOption("warn"); options(warn = -1)

  x <- try(as.numeric(x))

  options(warn = warn)

  if (any(is.na(x))) {
    FALSE
  }
  else {
    is.numeric(x) && all(abs(x - round(x)) < tol) 
  }
} ## is.wholenumber

ddirac <- function(x, location)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(x)) {
    if (all.equal(x[i], location) == TRUE) {
      f[i] = 1.0
    }
    else {
      f[i] = 0.0
    }
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## ddirac

RNGMIX <- function(Dataset = NULL,
  rseed = -1,
  n = NULL,
  Theta = NULL)
{
  message("RNGMIX Version 2.3.0");
  flush.console()

  if (is.null(Dataset)) {
    stop(sQuote("Dataset"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Dataset)) {
    stop(sQuote("Dataset"), " character vector is requested!", call. = FALSE)
  }

  if (!is.wholenumber(rseed)) {
    stop(sQuote("rseed"), " integer is requested!", call. = FALSE)
  }

  if (rseed > -1) {
    stop(sQuote("rseed"), " must be less than 0!", call. = FALSE)
  }

  if (is.null(n)) {
    stop(sQuote("n"), " must not be NULL!", call. = FALSE)
  }

  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer vector is requested!", call. = FALSE)
  }

  if (!all(n > 0)) {
    stop("all ", sQuote("n"), " must be greater than 0!", call. = FALSE)
  }

  if (is.null(Theta)) {
    stop(sQuote("Theta"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Theta)) {
    stop(sQuote("Theta"), " character matrix is requested!", call. = FALSE)
  }

  if (length(n) != ncol(Theta)) {
    stop("number of columns in ", sQuote("n"), " and ", sQuote("Theta"), " must match!", call. = FALSE)
  }

  C <- c("normal", "lognormal", "Weibull", "binomial", "Poisson", "Dirac")

  nrow <- nrow(Theta)
  ncol <- ncol(Theta)

  c <- as.integer(ncol)

  n <- rbind(n)

  Theta <- rbind(Theta)

  pdf <- array(data = NA, dim = c(nrow, ncol), dimnames = NULL)
  Theta1 <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)
  Theta2 <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)

  for (i in 1:ncol) {
    j <- 1; k <- 1
    
    while (j < nrow) {
      pdf[k, i] <- match.arg(Theta[j, i], C)

      if (pdf[k, i] %in% c("normal", "lognormal", "Weibull", "binomial")) {
        Theta1[k, i] <- as.numeric(Theta[j + 1, i])
        Theta2[k, i] <- as.numeric(Theta[j + 2, i])

        j <- j + 3; k <- k + 1
      }
      else
      if (pdf[k, i] %in% c("Poisson", "Dirac")) {
        Theta1[k, i] <- as.numeric(Theta[j + 1, i])
        Theta2[k, i] <- as.numeric(0.0)

        j <- j + 2; k <- k + 1
      }
    }
  }
  
  d <- k - 1; c <- ncol

  pdf <- pdf[1:d, ]; dim(pdf) <- c(d, c)
  Theta1 <- Theta1[1:d, ]; dim(Theta1) <- c(d, c)
  Theta2 <- Theta2[1:d, ]; dim(Theta2) <- c(d, c)

  xmin <- rbind(rep(+Inf, d))
  xmax <- rbind(rep(-Inf, d))

  IDum <- rseed

  RNGMIX <- list()

  RNGMIX$Dataset <- list()
    
  for (i in 1:length(Dataset)) {
    message("Dataset = ", Dataset[i])

    flush.console()

    output <- .C("RRNGMIX",
      IDum = as.integer(IDum), 
      d = as.integer(d),
      c = as.integer(c),
      N = as.integer(n),
      ParFamType = as.character(pdf),
      Par0 = as.double(Theta1),
      Par1 = as.double(Theta2),
      n = integer(1),
      X = double(sum(n) * d),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RNGMIX!", call. = FALSE); return(NA)
    }

    dim(output$X) <- c(output$n, d)

    xmin <- as.numeric(apply(rbind(xmin, output$X), 2, min))
    xmax <- as.numeric(apply(rbind(xmax, output$X), 2, max))

    RNGMIX$Dataset[[i]] <- as.data.frame(output$X, stringsAsFactors = FALSE)

    IDum <- IDum - 1
  }
  
  names(RNGMIX$Dataset) <- Dataset
  
  RNGMIX$w <- as.data.frame(rbind(format(n / output$n)), stringsAsFactors = FALSE)

  rownames(RNGMIX$w) <- "w"
  colnames(RNGMIX$w) <- paste("comp", if(c > 1) 1:c else "", sep = "")

  RNGMIX$Theta <- rbind(pdf, format(Theta1), format(Theta2))

  dim(RNGMIX$Theta) <- c(3 * d, c)

  if (d > 1) {
    rownames(RNGMIX$Theta) <- c(paste("pdf", 1:d, sep = ""),
      paste("theta1.", 1:d, sep = ""),
      paste("theta2.", 1:d, sep = ""))
  }
  else {
    rownames(RNGMIX$Theta) <- c("pdf", "theta1", "theta2")
  }

  Index <- NULL

  for (i in 1:d){
    Index <- c(Index, seq(from = i, to = i + 2 * d, by = d))
  }

  RNGMIX$Theta <- cbind(RNGMIX$Theta[Index, ])

  M <- match(RNGMIX$Theta[, 1], C)

  Index <- NULL

  for (i in 1:length(M)) {
    if (M[i] %in% c(5, 6)) {
      Index <- c(Index, i + 2)
    }
  }

  if (is.null(Index)) {
    RNGMIX$Theta <- as.data.frame(RNGMIX$Theta, stringsAsFactors = FALSE)
  }
  else {
    RNGMIX$Theta <- as.data.frame(RNGMIX$Theta[-Index, ], stringsAsFactors = FALSE)
  }

  colnames(RNGMIX$Theta) <- paste("comp", if(c > 1) 1:c else "", sep = "")

  RNGMIX$Variables <- rep("continuous", d)

  M <- na.omit(M)

  for (i in 1:length(M)) {
    if (M[i] %in% c(4, 5, 6)) {
      RNGMIX$Variables[i] <- "discrete"
    }
  }

  RNGMIX$ymin <- as.vector(xmin)
  RNGMIX$ymax <- as.vector(xmax)

  rm(list = ls()[!(ls() %in% c("RNGMIX"))])

  class(RNGMIX) <- "RNGMIX"
 
  return(RNGMIX)
} ## RNGMIX

print.RNGMIX <- function(x, ...) 
{
  if (missing(x) || (class(x) != "RNGMIX")) {
    stop(sQuote("x"), " object of class RNGMIX is requested!", call. = FALSE)
  }
  
  cat(paste("$w", "\n", sep = ""))
  
  print(x$w, ...) 
  
  cat(paste("\n", sep = ""))  
  
  cat(paste("$Theta", "\n", sep = ""))
  
  print(x$Theta, ...)  
  
  cat(paste("\n", sep = ""))  
  
  cat(paste("$ymin", "\n", sep = ""))
  
  print(x$ymin, ...)    

  cat(paste("\n", sep = ""))  
  
  cat(paste("$ymax", "\n", sep = ""))
  
  print(x$ymax, ...)    

  cat(paste("\n", sep = ""))  
  
  cat(paste("attr(,\"class\")", "\n", sep = ""))
  
  print(attr(x, "class"), ...)   

  invisible(x) 
} ## print.RNGMIX

.REBMIX <- function(Dataset = NULL, 
  Preprocessing = NULL, 
  D = 0.025, 
  cmax = 15,
  Criterion = "AIC",
  Variables = NULL,
  pdf = NULL,
  Theta1 = NULL,
  Theta2 = NULL,
  K = NULL,
  ymin = NULL,
  ymax = NULL,
  b = 1.0,
  ar = 0.1,
  Restraints = "loose")
{
  C <- c("normal", "lognormal", "Weibull", "binomial", "Poisson", "Dirac")

  REBMIX <- NULL
  REBMIX$Dataset <- Dataset
  REBMIX$w <- list()
  REBMIX$Theta <- list()
  REBMIX$Variables <- Variables
  REBMIX$summary <- list()
  REBMIX$pos <- 1

  for (i in 1:length(Dataset)) {
    DatasetName <- names(Dataset)[i]

    X <- as.matrix(Dataset[[i]])

    message("Dataset = ", DatasetName)
    flush.console()

    n <- nrow(X)
    d <- ncol(X)

    if (d < 1) {
      stop(sQuote("d"), " must be greater than 0!", call. = FALSE)
    }

    if (n < 1) {
      stop(sQuote("n"), " must be greater than 0!", call. = FALSE)
    }

    if (!is.null(Theta1)) {
      Theta1[is.na(Theta1)] <- 0
    }

    if (!is.null(Theta2)) {
      Theta2[is.na(Theta2)] <- 0
    }

    if (!is.null(ymin)) {
      ymin[is.na(ymin)] <- 0
    }

    if (!is.null(ymax)) {
      ymax[is.na(ymax)] <- 0
    }

    output <- .C("RREBMIX",
      PreType = as.character(Preprocessing), 
      D = as.double(D),
      cmax = as.integer(cmax),
      ICType = as.character(Criterion),
      d = as.integer(d),
      VarType = as.character(Variables),
      IniFamType = as.character(pdf),
      Ini0 = as.double(Theta1),
      Ini1 = as.double(Theta2),
      kmax = as.integer(length(K)),
      K = as.integer(K),
      ymin = as.double(ymin),
      ymax = as.double(ymax),
      b = as.double(b),
      ar = as.double(ar),
      ResType = as.character(Restraints),
      n = as.integer(n),
      X = as.double(X),
      k = integer(1),
      h = double(d),
      y0 = double(d),
      IC = double(1),
      logL = double(1),
      c = integer(1),
      W = double(cmax),        
      ParFamType = as.character(rep("THE_LONGEST_PARAMETRIC_FAMILY_TYPE", cmax * d)),
      Par0 = double(cmax * d),        
      Par1 = double(cmax * d),      
      ClcTime = integer(1),          
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in REBMIX!", call. = FALSE); return(NA)
    }

    c <- output$c

    length(output$h) <- d
    length(output$y0) <- d
    length(output$W) <- c
    length(output$ParFamType) <- c * d
    length(output$Par0) <- c * d
    length(output$Par1) <- c * d

    dim(output$ParFamType) <- c(d, c)
    dim(output$Par0) <- c(d, c)
    dim(output$Par1) <- c(d, c)

    REBMIX$w[[i]] <- as.data.frame(rbind(format(output$W)), stringsAsFactors = FALSE)

    rownames(REBMIX$w[[i]]) <- "w"
    colnames(REBMIX$w[[i]]) <- paste("comp", if(c > 1) 1:c else "", sep = "")

    REBMIX$Theta[[i]] <- rbind(output$ParFamType, format(output$Par0), format(output$Par1))

    dim(REBMIX$Theta[[i]]) <- c(3 * d, c)

    if (d > 1) {
      rownames(REBMIX$Theta[[i]]) <- c(paste("pdf", 1:d, sep = ""),
        paste("theta1.", 1:d, sep = ""),
        paste("theta2.", 1:d, sep = ""))
    }
    else {
      rownames(REBMIX$Theta[[i]]) <- c("pdf", "theta1", "theta2")
    }

    Index <- NULL

    for (j in 1:d){
      Index <- c(Index, seq(from = j, to = j + 2 * d, by = d))
    }

    REBMIX$Theta[[i]] <- cbind(REBMIX$Theta[[i]][Index, ])

    M <- match(REBMIX$Theta[[i]][, 1], C)

    Index <- NULL

    for (j in 1:length(M)) {
      if (M[j] %in% c(5, 6)) {
        Index <- c(Index, j + 2)
      }
    }

    if (is.null(Index)) {
      REBMIX$Theta[[i]] <- as.data.frame(REBMIX$Theta[[i]], stringsAsFactors = FALSE)
    }
    else {
      REBMIX$Theta[[i]] <- as.data.frame(REBMIX$Theta[[i]][-Index, ], stringsAsFactors = FALSE)
    }

    colnames(REBMIX$Theta[[i]]) <- paste("comp", if(c > 1) 1:c else "", sep = "")

    if (Preprocessing == "histogram") {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$PreType, 
        format(output$D),
        output$cmax, 
        output$ICType, 
        format(output$ar),
        output$ResType,
        output$c,
        format(output$b),
        output$k,
        format(output$y0),
        format(output$h),
        format(output$IC),
        format(output$logL))
    }
    else
    if (Preprocessing == "Parzen window") {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$PreType, 
        format(output$D),
        output$cmax, 
        output$ICType, 
        format(output$ar),
        output$ResType,
        output$c,
        format(output$b),
        output$k,
        format(output$h),
        format(output$IC),
        format(output$logL))
    }
    if (Preprocessing == "k-nearest neighbour") {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$PreType, 
        format(output$D),
        output$cmax, 
        output$ICType, 
        format(output$ar),
        output$ResType,
        output$c,
        format(output$b),
        output$k,
        format(output$h),
        format(output$IC),
        format(output$logL))
    }
  }

  REBMIX$summary <- as.data.frame(do.call("rbind", REBMIX$summary), stringsAsFactors = FALSE)

  if (Preprocessing == "histogram") {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "D", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "b",
      "v/k", 
      paste("y0", if (d > 1) 1:d else "", sep = ""), 
      paste("h", if (d > 1) 1:d else "", sep = ""), 
      "IC", 
      "logL")
  }
  else
  if (Preprocessing == "Parzen window") {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "D", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "b", 
      "v/k", 
      paste("h", if (d > 1) 1:d else "", sep = ""), 
      "IC", 
      "logL")
  }
  if (Preprocessing == "k-nearest neighbour") {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "D", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "b",
      "v/k", 
      paste("h", if (d > 1) 1:d else "", sep = ""),
      "IC", 
      "logL")
  }

  rm(list = ls()[!(ls() %in% c("REBMIX"))])

  class(REBMIX) <- "REBMIX"
 
  return(REBMIX)
} ## .REBMIX 

REBMIX <- function(Dataset = NULL, 
  Preprocessing = NULL, 
  D = 0.025, 
  cmax = 15,
  Criterion = "AIC",
  Variables = NULL,
  pdf = NULL,
  Theta1 = NULL,
  Theta2 = NULL,
  K = NULL,
  ymin = NULL,
  ymax = NULL,
  b = 1.0,
  ar = 0.1,
  Restraints = "loose")
{
  message("REBMIX Version 2.3.0");
  flush.console()

  if(is.null(Dataset)) {
    stop(sQuote("Dataset"), " must not be NULL!", call. = FALSE)
  }

  if (!is.list(Dataset)) {
    stop(sQuote("Dataset"), " list of data frames is requested!", call. = FALSE)
  }
  
  if (is.null(names(Dataset))) {
    names(Dataset) <- paste("dataset", 1:length(Dataset), sep = "")
  }
    
  for (i in 1:length(Dataset)) {
    if (!is.data.frame(Dataset[[i]])) {
      stop(sQuote("Dataset"), " list of data frames or character vector is requested!", call. = FALSE)
    }
      
    if ((is.na(names(Dataset)[i])) || (names(Dataset)[i] == "")) {
      names(Dataset)[i] <- paste("dataset", i, sep = "")  
    }
  }
  
  if(is.null(Preprocessing)) {
    stop(sQuote("Preprocessing"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Preprocessing)) {
    stop(sQuote("Preprocessing"), " character vector is requested!", call. = FALSE)
  }

  C <- c("histogram", "Parzen window", "k-nearest neighbour")

  Preprocessing <- match.arg(Preprocessing, C, several.ok = TRUE)

  if (!is.numeric(D)) {
    stop(sQuote("D"), " numeric is requested!", call. = FALSE)
  }

  if ((D < 0.0) || (D > 1.0)) {
    stop(sQuote("D"), " must be greater or equal than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

  if (!is.wholenumber(cmax)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (cmax < 1) {
    stop(sQuote("cmax"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character vector is requested!", call. = FALSE)
  }

  C <- c("AIC", "AIC3", "AIC4", "AICc", "AIC3", "BIC", "CAIC", "HQC", 
    "MDL2", "MDL5", "AWE", "CLC", "ICL", "ICL-BIC", "D", "SSE")

  Criterion <- match.arg(Criterion, C, several.ok = TRUE)

  if(is.null(Variables)) {
    stop(sQuote("Variables"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Variables)) {
    stop(sQuote("Variables"), " character vector is requested!", call. = FALSE)
  }

  C <- c("continuous", "discrete")

  Variables <- match.arg(Variables, C, several.ok = TRUE)

  if(is.null(pdf)) {
    stop(sQuote("pdf"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(pdf)) {
    stop(sQuote("pdf"), " character vector is requested!", call. = FALSE)
  }

  C <- c("normal", "lognormal", "Weibull", "binomial", "Poisson", "Dirac")

  pdf <- match.arg(pdf, C, several.ok = TRUE)

  if(is.null(K)) {
    stop(sQuote("K"), " must not be NULL!", call. = FALSE)
  }

  if(is.list(K)) {
    for (i in 1:length(K)) {
      if (!is.wholenumber(K[[i]])) {
        stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
      }

      if (!all(K[[i]] > 0)) {
        stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
      }
    }

    if (length(K) != length(Preprocessing)) {
      stop("lengths of ", sQuote("Preprocessing"), " and ", sQuote("K"), " must match!", call. = FALSE)
    }
  }
  else {
    if (!is.wholenumber(K)) {
      stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
    }

    if (!all(K > 0)) {
      stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
    }
  }

  if (!is.numeric(b)) {
    stop(sQuote("b"), " numeric is requested!", call. = FALSE)
  }

  if ((b < 0.0) || (b > 1.0)) {
    stop(sQuote("b"), " must be greater or equal than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

  if (!is.numeric(ar)) {
    stop(sQuote("ar"), " numeric is requested!", call. = FALSE)
  }

  if ((ar <= 0.0) || (ar > 1.0)) {
    stop(sQuote("ar"), " must be greater than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

  if (!is.character(Restraints)) {
    stop(sQuote("Restraints"), " character is requested!", call. = FALSE)
  }

  C <- c("rigid", "loose")

  Restraints <- match.arg(Restraints, C, several.ok = FALSE)

  REBMIX <- NULL
  REBMIX$Dataset <- Dataset
  REBMIX$w <- list()
  REBMIX$Theta <- list()
  REBMIX$Variables <- Variables
  REBMIX$summary <- NULL
  REBMIX$pos <- 1

  for (i in 1:length(Preprocessing)) {
    for (j in 1:length(Criterion)) {
      output <- .REBMIX(Dataset = Dataset, 
        Preprocessing = Preprocessing[i], 
        D = D, 
        cmax = cmax,
        Criterion = Criterion[j],
        Variables = Variables,
        pdf = pdf,
        Theta1 = Theta1,
        Theta2 = Theta2,
        K = if (is.list(K)) K[[i]] else K,
        ymin = ymin,
        ymax = ymax,
        b = b,
        ar = ar,
        Restraints = Restraints)

      for (k in (1:length(Dataset))) {
        REBMIX$w[[length(REBMIX$w) + 1]] <- output$w[[k]] 
        REBMIX$Theta[[length(REBMIX$Theta) + 1]] <- output$Theta[[k]]
      }

      if (is.null(REBMIX$summary)) {
        REBMIX$summary <- output$summary
      }
      else {
        REBMIX$summary <- merge(REBMIX$summary, output$summary, all = TRUE, sort = FALSE)
      }
    }
  }
  
  REBMIX$pos <- which(as.numeric(REBMIX$summary[, "logL"]) == max(as.numeric(REBMIX$summary[, "logL"])))  

  rm(list = ls()[!(ls() %in% c("REBMIX"))])

  class(REBMIX) <- "REBMIX"
 
  return(REBMIX)
} ## REBMIX

print.REBMIX <- function(x, pos = 1, ...) 
{
  if (missing(x) || (class(x) != "REBMIX")) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  if ((pos < 1) || (pos > nrow(x$summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x$summary), "!", call. = FALSE)
  }

  cat(paste("$w", "\n", sep = ""))
  
  print(x$w[[pos]], ...) 
  
  cat(paste("\n", sep = ""))  
  
  cat(paste("$Theta", "\n", sep = ""))
  
  print(x$Theta[[pos]], ...)  
  
  cat(paste("\n", sep = ""))  
  
  cat(paste("$summary", "\n", sep = ""))
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL"), names(x$summary), nomatch = 0)
  
  print(x$summary[pos, p], ...)  

  cat(paste("\n", sep = ""))  
  
  cat(paste("attr(,\"class\")", "\n", sep = ""))
  
  print(attr(x, "class"), ...)   

  invisible(x) 
} ## print.REBMIX 

coef.REBMIX <- function(object, pos = 1, ...) 
{
  if (missing(object) || (class(object) != "REBMIX")) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  print(rbind(object$w[[pos]], object$Theta[[pos]]), ...)
} ## coef.REBMIX

summary.REBMIX <- function(object, ...) 
{
  if (missing(object) || (class(object) != "REBMIX")) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL"), names(object$summary), nomatch = 0)
  
  print(object$summary[p], ...)
  
  cat(paste("Maximum logL = ", object$summary[object$pos, "logL"], " at pos = ", object$pos, ".\n", sep = "", collapse = "")) 
} ## summary.REBMIX

.densKNearestNeighbour.x <- function(x, k, hx) 
{
  output <- .C("RdensKNearestNeighbourX",
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    k = as.integer(k),
    hx = as.double(hx),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densKNearestNeighbour.x!", call. = FALSE); return(NA)
  }

  i <- order(output$y)

  output$x <- output$x[i] 
  output$y <- output$y[i] 

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densKNearestNeighbour.x

.densKNearestNeighbour.xy <- function(x, y, k, hx, hy) 
{
  output <- .C("RdensKNearestNeighbourXY",
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    k = as.integer(k),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densKNearestNeighbour.xy!", call. = FALSE); return(NA)
  }

  i <- order(output$z)

  output$x <- output$x[i] 
  output$y <- output$y[i] 
  output$z <- output$z[i] 

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
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densParzenWindow.x!", call. = FALSE); return(NA)
  }

  i <- order(output$y)

  output$x <- output$x[i] 
  output$y <- output$y[i] 
  
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
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densParzenWindow.xy!", call. = FALSE); return(NA)
  }

  i <- order(output$z)

  output$x <- output$x[i] 
  output$y <- output$y[i] 
  output$z <- output$z[i] 

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densParzenWindow.xy

.densHistogram.x <- function(k, x, x0, hx, cx) 
{
  output <- .C("RdensHistogramX",
    k = as.integer(k),
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    x0 = as.double(x0),
    hx = as.double(hx),
    cx = as.integer(toupper(cx) == "DISCRETE"),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densHistogram.x!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.x 

.densHistogram.xy <- function(k, x, y, x0, y0, hx, hy, cx, cy) 
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
    cx = as.integer(toupper(cx) == "DISCRETE"),
    cy = as.integer(toupper(cy) == "DISCRETE"),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densHistogram.xy!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k
  length(output$z) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.xy

.dfmix.x <- function(x, w, xTheta, ...) 
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    C <- toupper(xTheta[[i]]$pdf)

    if (C == "NORMAL") {
      fix <- dnorm(as.numeric(x), mean = as.numeric(xTheta[[i]]$theta1), sd = as.numeric(xTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "LOGNORMAL") {
      fix <- dlnorm(as.numeric(x), meanlog = as.numeric(xTheta[[i]]$theta1), sdlog = as.numeric(xTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "WEIBULL") {
      fix <- dweibull(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "BINOMIAL") {
      fix <- dbinom(as.integer(x), size = as.integer(xTheta[[i]]$theta1), prob = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else 
    if (C == "POISSON") {
      fix <- dpois(as.integer(x), lambda = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else 
    if (C == "DIRAC") {
      fix <- ddirac(x, location = as.numeric(xTheta[[i]]$theta1))
    }

    f <- f + w[i] * fix
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## .dfmix.x

.dfmix.xy <- function(x, y, w, xTheta, yTheta, ...) 
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    C <- toupper(xTheta[[i]]$pdf)

    if (C == "NORMAL") {
      fix <- dnorm(as.numeric(x), mean = as.numeric(xTheta[[i]]$theta1), sd = as.numeric(xTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "LOGNORMAL") {
      fix <- dlnorm(as.numeric(x), meanlog = as.numeric(xTheta[[i]]$theta1), sdlog = as.numeric(xTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "WEIBULL") {
      fix <- dweibull(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "BINOMIAL") {
      fix <- dbinom(as.integer(x), size = as.integer(xTheta[[i]]$theta1), prob = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else 
    if (C == "POISSON") {
      fix <- dpois(as.integer(x), lambda = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else
    if (C == "DIRAC") {
      fix <- ddirac(x, location = as.numeric(xTheta[[i]]$theta1))
    }

    C <- toupper(yTheta[[i]]$pdf)

    if (C == "NORMAL") {
      fiy <- dnorm(as.numeric(y), mean = as.numeric(yTheta[[i]]$theta1), sd = as.numeric(yTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "LOGNORMAL") {
      fiy <- dlnorm(as.numeric(y), meanlog = as.numeric(yTheta[[i]]$theta1), sdlog = as.numeric(yTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "WEIBULL") {
      fiy <- dweibull(as.numeric(y), scale = as.numeric(yTheta[[i]]$theta1), shape = as.numeric(yTheta[[i]]$theta2), ...)
    } 
    else 
    if (C == "BINOMIAL") {
      fiy <- dbinom(as.integer(y), size = as.integer(yTheta[[i]]$theta1), prob = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else 
    if (C == "POISSON") {
      fiy <- dpois(as.integer(y), lambda = as.numeric(yTheta[[i]]$theta1), ...)
    }
    else
    if (C == "DIRAC") {
      fiy <- ddirac(y, location = as.numeric(yTheta[[i]]$theta1))
    }

    f <- f + w[i] * fix * fiy
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## .dfmix.xy

.extractTheta <- function(Theta) 
{
  C <- c("NORMAL", "LOGNORMAL", "WEIBULL", "BINOMIAL", "POISSON", "DIRAC")

  nrow <- nrow(Theta)
  ncol <- ncol(Theta)

  output <- array(data = list(NULL), dim = c(nrow, ncol), dimnames = NULL)

  for (i in 1:ncol) {
    M <- match(toupper(Theta[, i]), C)

    j <- 1;

    for (k in 1:length(M)) {
      if (M[k] %in% c(1, 2, 3, 4)) {
        output[[j, i]]$pdf <- Theta[k, i]
        output[[j, i]]$theta1 <- as.numeric(Theta[k + 1, i])
        output[[j, i]]$theta2 <- as.numeric(Theta[k + 2, i])

        j <- j + 1
      }
      else
      if (M[k] %in% c(5, 6)) {
        output[[j, i]]$pdf <- Theta[k, i]
        output[[j, i]]$theta1 <- as.numeric(Theta[k + 1, i])

        j <- j + 1
      }
    }
  }

  j <- j - 1

  output <- output[1:j, ]; dim(output) <- c(j, ncol)

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .extractTheta

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
  contour.drawlabels = FALSE,
  contour.labcex = 0.8, 
  contour.method = "flattest",
  contour.nlevels = 12, ...) 
{
  if (missing(x) || (class(x) != "REBMIX")) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  if ((pos < 1) || (pos > nrow(x$summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x$summary), "!", call. = FALSE)
  }

  if (!is.wholenumber(nrow)) {
    stop(sQuote("nrow"), " integer is requested!", call. = FALSE)
  }

  if(nrow < 1) {
    stop(sQuote("nrow"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(ncol)) {
    stop(sQuote("ncol"), " integer is requested!", call. = FALSE)
  }

  if(ncol < 1) {
    stop(sQuote("ncol"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(npts)) {
    stop(sQuote("npts"), " integer is requested!", call. = FALSE)
  }

  if(npts < 1) {
    stop(sQuote("npts"), " must be greater than 0!", call. = FALSE)
  }

  ni <- ncol(x$summary)

  Theta <- .extractTheta(x$Theta[[pos]])

  d <- nrow(Theta)

  if (d == 1) {
    nrow <- 1; ncol <- 1
  }
  else {
    nrow <- max(1, nrow)
    ncol <- max(1, ncol)
  }
    
  n <- d * (d - 1) / 2

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
    
  C <- x$summary[pos, "Preprocessing"]
  
  if (.Device == "tikz output") {
    item <- list()

    item[[1]] <- "$\\textrm{Dataset}$"
    item[[2]] <- "$\\; = \\;$"
    item[[3]] <- paste("$\\textrm{", x$summary[pos, "Dataset"], "}$, ", sep = "")

    item[[4]] <- "$\\textrm{Preprocessing}$"
    item[[5]] <- "$\\; = \\;$"
    item[[6]] <- paste("$\\textrm{", x$summary[pos, "Preprocessing"], "}$, ", sep = "")

    item[[7]] <- "$\\textrm{Restraints}$"
    item[[8]] <- "$\\; = \\;$"
    item[[9]] <- paste("$\\textrm{", x$summary[pos, "Restraints"], "}$, ", sep = "")

    item[[10]] <- "$D$"
    item[[11]] <- "$\\; = \\;$"
    item[[12]] <- paste("$", x$summary[pos, "D"], "$, ", sep = "")

    item[[13]] <- "$c_{\\mathrm{max}}$"
    item[[14]] <- "$\\; = \\;$"
    item[[15]] <- paste("$", x$summary[pos, "cmax"], "$, ", sep = "")

    item[[16]] <- "$a_{\\mathrm{r}}$"
    item[[17]] <- "$\\; = \\;$"
    item[[18]] <- paste("$", x$summary[pos, "ar"], "$, ", sep = "")

    item[[19]] <- "$c$"
    item[[20]] <- "$\\; = \\;$"
    item[[21]] <- paste("$", x$summary[pos, "c"], "$, ", sep = "")

    item[[22]] <- "$b$"
    item[[23]] <- "$\\; = \\;$"
    item[[24]] <- paste("$", x$summary[pos, "b"], "$, ", sep = "")

    if (C == "histogram") {
      item[[25]] <- "$v$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x$summary[pos, "v/k"], "$, ", sep = "")    
    } 
    else 
    if (C == "Parzen window") {
      item[[25]] <- "$v$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x$summary[pos, "v/k"], "$, ", sep = "")
    } 
    else
    if (C == "k-nearest neighbour") {
      item[[25]] <- "$k$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x$summary[pos, "v/k"], "$, ", sep = "")
    }

    item[[28]] <- paste("$\\mathrm{", x$summary[pos, "Criterion"], "}$", sep = "")
    item[[29]] <- "$\\; = \\;$"
    item[[30]] <- paste("$", format(x$summary[pos, "IC"]), "$, ", sep = "")

    item[[31]] <- "$\\mathrm{log}\\, L$"
    item[[32]] <- "$\\; = \\;$"
    item[[33]] <- paste("$", format(x$summary[pos, "logL"]), "$.", sep = "")

    i <- 1; legend <- list(); legend[[i]] <- item[[1]]

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
    item <- list()

    item[[1]] <- "Dataset"
    item[[2]] <- " = "
    item[[3]] <- paste(x$summary[pos, "Dataset"], ", ", sep = "")

    item[[4]] <- "Preprocessing"
    item[[5]] <- " = "
    item[[6]] <- paste(x$summary[pos, "Preprocessing"], ", ", sep = "")

    item[[7]] <- "Restraints"
    item[[8]] <- " = "
    item[[9]] <- paste(x$summary[pos, "Restraints"], ", ", sep = "")

    item[[10]] <- "D"
    item[[11]] <- " = "
    item[[12]] <- paste(x$summary[pos, "D"], ", ", sep = "")

    item[[13]] <- bquote(c[max])
    item[[14]] <- " = "
    item[[15]] <- paste(x$summary[pos, "cmax"], ", ", sep = "")

    item[[16]] <- bquote(a[r])
    item[[17]] <- " = "
    item[[18]] <- paste(x$summary[pos, "ar"], ", ", sep = "")

    item[[19]] <- "c"
    item[[20]] <- " = "
    item[[21]] <- paste(x$summary[pos, "c"], ", ", sep = "")

    item[[22]] <- "b"
    item[[23]] <- " = "
    item[[24]] <- paste(x$summary[pos, "b"], ", ", sep = "")

    if (C == "histogram") {
      item[[25]] <- "v"
      item[[26]] <- " = "
      item[[27]] <- paste(x$summary[pos, "v/k"], ", ", sep = "")
    } 
    else 
    if (C == "Parzen window") {
      item[[25]] <- "v"
      item[[26]] <- " = "
      item[[27]] <- paste(x$summary[pos, "v/k"], ", ", sep = "")
    } 
    else
    if (C == "k-nearest neighbour") {
      item[[25]] <- "k"
      item[[26]] <- " = "
      item[[27]] <- paste(x$summary[pos, "v/k"], ", ", sep = "")
    }

    item[[28]] <- as.character(x$summary[pos, "Criterion"])
    item[[29]] <- " = "
    item[[30]] <- paste(format(x$summary[pos, "IC"]), ", ", sep = "")

    item[[31]] <- "log L"
    item[[32]] <- " = "
    item[[33]] <- paste(format(x$summary[pos, "logL"]), ".", sep = "")

    i <- 1; legend <- list(); legend[[i]] <- bquote(.(item[[1]]))

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

  Dataset <- as.character(x$summary[pos, "Dataset"])

  ey <- as.matrix(x$Dataset[[which(names(x$Dataset) == x$summary[pos, "Dataset"])]])

  y0 <- array(data = 0.0, dim = d, dimnames = NULL)
  h <- array(data = 0.0, dim = d, dimnames = NULL)

  lim <- array(data = 0.0, dim = c(2, d), dimnames = NULL)

  py <- list(d)

  Variables <- x$Variables

  b <- as.numeric(x$summary[pos, "b"])

  for (i in 1:d) {
    if (C == "histogram") {
      k <- as.numeric(x$summary[pos, "v/k"])
      y0[i] <- as.numeric(x$summary[pos, paste("y0", if (d > 1) i, sep = "")])
      h[i] <- as.numeric(x$summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    } 
    else 
    if (C == "Parzen window") {
      h[i] <- as.numeric(x$summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    } 
    else
    if (C == "k-nearest neighbour") {
      k <- as.numeric(x$summary[pos, "v/k"])

      h[i] <- as.numeric(x$summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }

    if (Variables[i] == "discrete") {
      py[[i]] <- seq(from = lim[1, i], to = lim[2, i], by = 1.0)
    }
    else {
      py[[i]] <- seq(from = lim[1, i], to = lim[2, i], length.out = npts / ncol)
    }
  }

  w <- as.numeric(x$w[[pos]])

  if (n > 0) {
    ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
      space = "rgb",
      interpolate = "linear")

    figno <- 0

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        if (C == "histogram") {
          edens <- .densHistogram.xy(k, ey[, i], ey[, j], y0[i], y0[j], h[i], h[j], Variables[i], Variables[j]) 
        } 
        else 
        if (C == "Parzen window") {
          edens <- .densParzenWindow.xy(ey[, i], ey[, j], h[i], h[j])
        } 
        else
        if (C == "k-nearest neighbour") {
          edens <- .densKNearestNeighbour.xy(ey[, i], ey[, j], k, h[i], h[j])
        }

        pdens <- outer(py[[i]], py[[j]], ".dfmix.xy", w, Theta[i, ], Theta[j, ])

        zlim <- range(edens$z, finite = TRUE); zmax <- max(zlim[2], pdens)

        zlim <- zlim / zmax
   
        plot(x = edens$x, 
          y = edens$y, 
          type = "p",
          main = "", 
          sub = "",   
          xlab = "",
          ylab = "",  
          col = rgb(ramp(edens$z / zmax), maxColorValue = 255),
          axes = FALSE,
          lwd = 1,
          cex = plot.cex,
          pch = plot.pch)

        if ((Variables[i] == "discrete") && (Variables[j] == "discrete")) {
          points(x = rep(py[[i]], length(py[[j]])),
            y = rep(py[[j]], each = length(py[[i]])), 
            type = "p", 
            xlab = "",
            ylab = "", 
            col = rgb(ramp(pdens / zmax), maxColorValue = 255),
            lwd = 1,
            cex = plot.cex * 0.5,
            pch = plot.pch)
        }
        else
        if ((Variables[i] == "discrete") && (Variables[j] == "continuous")) {
          for (l in 1:length(py[[i]])) {
            tx <- rep(py[[i]][l], length(py[[j]]))
            ty <- py[[j]]

            s <- 1:(length(tx)-1)
            
            segments(x0 = tx[s], 
              y0 = ty[s], 
              x1 = tx[s + 1], 
              y1 = ty[s + 1], 
              xlab = "",
              ylab = "", 
              col = rgb(ramp((pdens[l, s] + pdens[l, s + 1]) / zmax / 2.0), maxColorValue = 255),
              cex = plot.cex)
          }
        }
        else 
        if ((Variables[i] == "continuous") && (Variables[j] == "discrete")) {
          for (l in 1:length(py[[j]])) {
            tx <- py[[i]]
            ty <- rep(py[[j]][l], length(py[[i]]))

            s <- 1:(length(tx)-1)
            
            segments(x0 = tx[s], 
              y0 = ty[s], 
              x1 = tx[s + 1], 
              y1 = ty[s + 1], 
              xlab = "",
              ylab = "", 
              col = rgb(ramp((pdens[s, l] + pdens[s + 1, l]) / zmax / 2.0), maxColorValue = 255),
              cex = plot.cex)
          }
        }
        else {
          levels <- 10^seq(from = log(zlim[1]), to = log(zlim[2]), length.out = contour.nlevels)

          contour(x = py[[i]],
            y = py[[j]],
            z = pdens / zmax, 
            levels = levels,
            xlim = lim[, i],
            ylim = lim[, j],
            zlim = zlim,
            labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
            axes = FALSE, frame.plot = FALSE,
            col = rgb(ramp(levels), maxColorValue = 255), 
            add = TRUE)
        }
   
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
          for (l in 1:length(legend)) {
            mtext(text = legend[[l]],
              side = 1, 
              line = l - 1, 
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
    if (C == "histogram") {
      edens <- .densHistogram.x(k, ey[, 1], y0[1], h[1], Variables[1])
    } 
    else 
    if (C == "Parzen window") {
      edens <- .densParzenWindow.x(ey[, 1], h[1])
    } 
    else
    if (C == "k-nearest neighbour") {
      edens <- .densKNearestNeighbour.x(ey[, 1], k, h[1])
    }

    pdens <- .dfmix.x(py[[1]], w, Theta[1, ])
        
    ylim <- c(0.0, max(edens$y, pdens))

    plot(x = edens$x, 
      y = edens$y, 
      type = "p",
      main = "", 
      sub = "",   
      xlab = "",
      ylab = "",  
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

    points(x = py[[1]], 
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

    for (l in 1:length(legend)) {
      mtext(text = legend[[l]],
        side = 1, 
        line = l - 1, 
        outer = TRUE,
        adj = 0.5,
        padj = 0.0,
        cex = cex)
    }
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])
  
  invisible(opar)
} ## plot.REBMIX

predict.list <- function(object,
  P = NULL,
  Dataset = NULL, ...) 
{
  if (missing(object)) {
    stop(sQuote("object"), " object or list of classes REBMIX is requested!", call. = FALSE)
  }

  if (class(object) == "list") {
    o <- length(object)

    for (i in 1:o) {
      if (class(object[[i]]) != "REBMIX") {
        stop(sQuote("object"), " list of classes REBMIX is requested!", call. = FALSE)
      }
    }
  }
  else
  if (class(object) == "REBMIX") {
    o <- 1;

    object[[1]] <- object
  }
  else {
    stop(sQuote("object"), " object or list of classes REBMIX is requested!", call. = FALSE)
  }

  if(is.null(Dataset)) {
    stop(sQuote("Dataset"), " must not be NULL!", call. = FALSE)
  }

  Dataset <- as.data.frame(Dataset)

  s <- nrow(object[[1]]$summary)

  if (o > 1) {
    for (i in 2:o) {
      if (s != nrow(object[[i]]$summary)) {
        stop(sQuote("object"), " list of classes REBMIX with equal number of classes is requested!", call. = FALSE)
      }
    }
  }

  if (s > 1) {
    if (is.null(P)) {
      stop(sQuote("P"), " must not be NULL!", call. = FALSE)
    }

    if (s != length(P)) {
      stop(sQuote("object"), " and ", sQuote("P"), " must be of the same length!", call. = FALSE)
    }
  }

  d <- array(data = NA, dim = o, dimnames = NULL)

  c <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  w <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  pdf <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  Theta1 <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  Theta2 <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  C <- c("NORMAL", "LOGNORMAL", "WEIBULL", "BINOMIAL", "POISSON", "DIRAC")

  for (io in 1:o) {
    for (is in 1:s) {
      nrow <- nrow(object[[io]]$Theta[[is]])
      ncol <- ncol(object[[io]]$Theta[[is]])

      c[[io, is]] <- as.integer(ncol)

      w[[io, is]] <- as.numeric(object[[io]]$w[[is]])

      pdf[[io, is]] <- array(data = NA, dim = c(nrow, ncol), dimnames = NULL)
      Theta1[[io, is]] <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)
      Theta2[[io, is]] <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)

      for (j in 1:ncol) {
        M <- match(toupper(object[[io]]$Theta[[is]][, j]), C)

        d[io] <- 1;

        for (l in 1:length(M)) {
          if (M[l] %in% c(1, 2, 3, 4)) {
            pdf[[io, is]][d[io], j] <- object[[io]]$Theta[[is]][l, j]
            Theta1[[io, is]][d[io], j] <- as.numeric(object[[io]]$Theta[[is]][l + 1, j])
            Theta2[[io, is]][d[io], j] <- as.numeric(object[[io]]$Theta[[is]][l + 2, j])

            d[io] <- d[io] + 1
          }
          else
          if (M[l] %in% c(5, 6)) {
            pdf[[io, is]][d[io], j] <- object[[io]]$Theta[[is]][l, j]
            Theta1[[io, is]][d[io], j] <- as.numeric(object[[io]]$Theta[[is]][l + 1, j])

            d[io] <- d[io] + 1
          }
        }
      }

      d[io] <- d[io] - 1

      pdf[[io, is]] <- pdf[[io, is]][1:d[io], ]; dim(pdf[[io, is]]) <- c(d[io], ncol)
      Theta1[[io, is]] <- Theta1[[io, is]][1:d[io], ]; dim(Theta1[[io, is]]) <- c(d[io], ncol)
      Theta2[[io, is]] <- Theta2[[io, is]][1:d[io], ]; dim(Theta2[[io, is]]) <- c(d[io], ncol)
    }
  }

  if (s > 1) {
    message("RCLSMIX Version 2.3.0");
    flush.console()

    output <- .C("RCLSMIX",
      n = as.integer(nrow(Dataset)),
      X = as.double(unlist(Dataset)),
      s = as.integer(s),
      o = as.integer(o),
      d = as.integer(d),
      c = as.integer(unlist(c)),
      W = as.double(unlist(w)),
      ParFamType = as.character(unlist(pdf)),
      Par0 = as.double(unlist(Theta1)),
      Par1 = as.double(unlist(Theta2)),
      P = as.double(unlist(P)),
      Z = integer(nrow(Dataset)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLSMIX!", call. = FALSE); return(NA)
    }
  }

  output <- as.factor(output$Z)
  levels(output) <- 0:(length(P) -1)

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## predict.list
