.pfmix.x <- function(x, w, xTheta, ...)
{
  n <- length(x)

  f <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:length(w)) {
    if (xTheta[[i]]$pdf == .rebmix$pdf[1]) {
      fix <- pnorm(as.numeric(x), mean = as.numeric(xTheta[[i]]$theta1), sd = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[2]) {
      fix <- plnorm(as.numeric(x), meanlog = as.numeric(xTheta[[i]]$theta1), sdlog = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[3]) {
      fix <- pweibull(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[4]) {
      fix <- pbinom(as.integer(x), size = as.integer(xTheta[[i]]$theta1), prob = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[5]) {
      fix <- ppois(as.integer(x), lambda = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[6]) {
      fix <- pdirac(as.numeric(x), location = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[7]) {
      fix <- pgamma(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[8]) {
      fix <- punif(as.numeric(x), min = as.numeric(xTheta[[i]]$theta1), max = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else    
    if (xTheta[[i]]$pdf == .rebmix$pdf[9]) {
      output <- .C(C_RvonMisesCdf,
        n = as.integer(n),
        y = as.double(x),
        Mean = as.double(xTheta[[i]]$theta1),
        Kappa = as.double(xTheta[[i]]$theta2),
        F = double(n),
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
         
      fix <- output$F
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[10]) {
      output <- .C(C_RGumbelCdf,
        n = as.integer(n),
        y = as.double(x),
        Mean = as.double(xTheta[[i]]$theta1),
        Sigma = as.double(xTheta[[i]]$theta2),
        Xi = as.double(xTheta[[i]]$theta3),
        F = double(n),
        PACKAGE = "rebmix")

      fix <- output$F    
    } 

    f <- f + w[i] * fix
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .pfmix.x
