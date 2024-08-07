setMethod(".IC",
          signature(x = "REBMIX"),
function(x, Criterion, pos, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character vector is requested!", call. = FALSE)
  }

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = TRUE)

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }

  Dataset <- as.character(x@summary[pos, "Dataset"])
  
  Dataset <- x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]
  
  if (as.character(class(Dataset)) == "data.frame") {
    Y.type <- 0
    
    X <- as.matrix(Dataset)

    n <- nrow(X)
    d <- ncol(X)    
  }
  else
  if (as.character(class(Dataset)) == "Histogram") {
    Y.type <- 1
    
    X <- as.matrix(Dataset@Y)

    n <- nrow(X)
    d <- ncol(X) - 1 
  }  

  h <- as.double(x@summary[pos, paste("h", if (d > 1) 1:d else "", sep = "")])

  Names <- names(x@Theta[[pos]])

  c <- length(x@w[[pos]])

  pdf <- unlist(x@Theta[[pos]][grep("pdf", Names)])

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

  theta1 <- unlist(x@Theta[[pos]][grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(x@Theta[[pos]][grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0
  
  theta3 <- unlist(x@Theta[[pos]][grep("theta3", Names)])

  theta3[is.na(theta3)] <- 0  

  length(pdf) <- d

  C <- x@summary[pos, "Preprocessing"]
  
  if (Y.type == 0) {
    if (is.na(C)) {
      output <- .C(C_RInformationCriterionMIX,
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        Theta1 = as.double(theta1),
        Theta2 = as.double(theta2),
        Theta3 = as.double(theta3),        
        length.pdf = as.integer(d),
        length.Theta = as.integer(3),
        length.theta = as.integer(c(d, d, d)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2, theta3)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
    
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }    
    }
    else
    if (C == .rebmix$Preprocessing[1]) {
      y0 <- as.double(x@summary[pos, paste("y0", if (d > 1) 1:d else "", sep = "")])
      ymin <- as.double(x@summary[pos, paste("ymin", if (d > 1) 1:d else "", sep = "")])
      ymax <- as.double(x@summary[pos, paste("ymax", if (d > 1) 1:d else "", sep = "")])

      output <- .C(C_RInformationCriterionHMIX,
        h = as.double(h),
        y0 = as.double(y0),
        ymin = as.double(ymin),
        ymax = as.double(ymax),
        k = as.integer(x@summary[pos, "v/k"]),
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        Theta1 = as.double(theta1),
        Theta2 = as.double(theta2),
        Theta3 = as.double(theta3),        
        length.pdf = as.integer(d),
        length.Theta = as.integer(3),
        length.theta = as.integer(c(d, d, d)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2, theta3)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
    
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }   
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      output <- .C(C_RInformationCriterionKDEMIX,
        h = as.double(h),
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        Theta1 = as.double(theta1),
        Theta2 = as.double(theta2),
        Theta3 = as.double(theta3),        
        length.pdf = as.integer(d),
        length.Theta = as.integer(3),
        length.theta = as.integer(c(d, d, d)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2, theta3)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
    
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }   
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      k <- as.integer(x@summary[pos, "v/k"])

      output <- .C(C_RInformationCriterionKNNMIX,
        h = as.double(h),
        k = as.integer(x@summary[pos, "v/k"]),
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        Theta1 = as.double(theta1),
        Theta2 = as.double(theta2),
        Theta3 = as.double(theta3),        
        length.pdf = as.integer(d),
        length.Theta = as.integer(3),
        length.theta = as.integer(c(d, d, d)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2, theta3)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
    
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }   
    }
  }  
  else
  if (Y.type == 1) {
    output <- .C(C_RInformationCriterionKMIX,
      h = as.double(h),
      Criterion = as.character(Criterion),
      c = as.integer(c),
      w = as.double(x@w[[pos]]),
      Theta1 = as.double(theta1),
      Theta2 = as.double(theta2),
      Theta3 = as.double(theta3),      
      length.pdf = as.integer(d),
      length.Theta = as.integer(3),
      length.theta = as.integer(c(d, d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2, theta3)),
      n = as.integer(n),
      x = as.double(X),
      IC = double(1),
      logL = double(1),
      M = integer(1),
      D = double(1),
      error = integer(9),
      PACKAGE = "rebmix")

    error <- error.to.string(output$error);
      
    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
   
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    }    
  }  

  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)
}) ## .IC

setMethod(".IC",
          signature(x = "REBMVNORM"),
function(x, Criterion, pos, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMVNORM is requested!", call. = FALSE)
  }

  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character vector is requested!", call. = FALSE)
  }

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = TRUE)

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }

  Dataset <- as.character(x@summary[pos, "Dataset"])
  
  Dataset <- x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]

  if (as.character(class(Dataset)) == "data.frame") {
    Y.type <- 0
    
    X <- as.matrix(Dataset)

    n <- nrow(X)
    d <- ncol(X)    
  }
  else
  if (as.character(class(Dataset)) == "Histogram") {
    Y.type <- 1
    
    X <- as.matrix(Dataset@Y)

    n <- nrow(X)
    d <- ncol(X) - 1 
  } 

  h <- as.double(x@summary[pos, paste("h", if (d > 1) 1:d else "", sep = "")])

  Names <- names(x@Theta[[pos]])

  c <- length(x@w[[pos]])

  pdf <- unlist(x@Theta[[pos]][grep("pdf", Names)])

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

  theta1 <- unlist(x@Theta[[pos]][grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(x@Theta[[pos]][grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0

  length(pdf) <- d

  C <- x@summary[pos, "Preprocessing"]
  
  if (Y.type == 0) {
    if (is.na(C)) {
      output <- .C(C_RInformationCriterionMVNORM,
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        length.pdf = as.integer(d),
        length.Theta = as.integer(4),
        length.theta = as.integer(c(d, d * d, -d * d, -1)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
    
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }       
    }
    else  
    if (C == .rebmix$Preprocessing[1]) {
      y0 <- as.double(x@summary[pos, paste("y0", if (d > 1) 1:d else "", sep = "")])
      ymin <- as.double(x@summary[pos, paste("ymin", if (d > 1) 1:d else "", sep = "")])
      ymax <- as.double(x@summary[pos, paste("ymax", if (d > 1) 1:d else "", sep = "")])

      output <- .C(C_RInformationCriterionHMVNORM,
        h = as.double(h),
        y0 = as.double(y0),
        ymin = as.double(ymin),
        ymax = as.double(ymax),
        k = as.integer(x@summary[pos, "v/k"]),
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        length.pdf = as.integer(d),
        length.Theta = as.integer(4),
        length.theta = as.integer(c(d, d * d, -d * d, -1)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
    
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }   
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      output <- .C(C_RInformationCriterionKDEMVNORM,
        h = as.double(h),
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        length.pdf = as.integer(d),
        length.Theta = as.integer(4),
        length.theta = as.integer(c(d, d * d, -d * d, -1)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

     error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
    
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }   
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      k <- as.integer(x@summary[pos, "v/k"])

      output <- .C(C_RInformationCriterionKNNMVNORM,
        h = as.double(h),
        k = as.integer(x@summary[pos, "v/k"]),
        Criterion = as.character(Criterion),
        c = as.integer(c),
        w = as.double(x@w[[pos]]),
        length.pdf = as.integer(d),
        length.Theta = as.integer(4),
        length.theta = as.integer(c(d, d * d, -d * d, -1)),
        pdf = as.character(pdf),
        Theta = as.double(c(theta1, theta2)),
        n = as.integer(n),
        x = as.double(X),
        IC = double(1),
        logL = double(1),
        M = integer(1),
        D = double(1),
        error = integer(9),
        PACKAGE = "rebmix")

      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
   
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      } 
    }
  }
  else
  if (Y.type == 1) {
    output <- .C(C_RInformationCriterionKMVNORM,
      h = as.double(h),
      Criterion = as.character(Criterion),
      c = as.integer(c),
      w = as.double(x@w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, -d * d, -1)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      IC = double(1),
      logL = double(1),
      M = integer(1),
      D = double(1),
      error = integer(9),
      PACKAGE = "rebmix")

    error <- error.to.string(output$error);
      
    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
   
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    } 
  }    

  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)
}) ## .IC

setMethod("logL",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC", pos = pos, ...)

  output <- output$logL

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## logL

setMethod("AIC",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## AIC

setMethod("AIC3",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC3", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## AIC3

setMethod("AIC4",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC4", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## AIC4

setMethod("AICc",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AICc", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## AICc

setMethod("BIC",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "BIC", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## BIC

setMethod("CAIC",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "CAIC", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## CAIC

setMethod("HQC",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "HQC", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## HQC

setMethod("MDL2",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "MDL2", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## MDL2

setMethod("MDL5",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "MDL5", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## MDL5

setMethod("AWE",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AWE", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## AWE

setMethod("CLC",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "CLC", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## CLC

setMethod("ICL",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "ICL", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## ICL

setMethod("ICLBIC",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "ICL-BIC", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## ICLBIC

setMethod("PRD",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "D", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## PRD

setMethod("SSE",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "SSE", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## SSE

setMethod("PC",
          signature(x = "ANY"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "PC", pos = pos, ...)

  output <- output$IC

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## PC
