
nonlinEffect <- function(fit, select=NULL, renames=NULL) 
{
  # Check if the model contains nonlinear variables.
  if(length(fit$nonlinear)==0){
    cat("There are no nonlinear effects in the model.  There will be no plots.\n")
    nonlineffs <- NA
    class(nonlineffs) <- "plaqreffect"
    return(nonlineffs)
  }

  plotvars <- fit$nonlinear
  # if select, only use selected nonlinear variables
  if(!is.null(select)){
    plotvars <- unique(plotvars[match(select, plotvars)])
    plotvars <- plotvars[!is.na(plotvars)]

    # Check if at least one nonlinear variable is selected, otherwise give feedback.
    if(all(is.na(plotvars))){
      cat("No nonlinear effects are selected.  There will be no plots.\n")
      nonlineffs <- NA
      class(nonlineffs) <- "plaqreffect"
      return(nonlineffs)
    }
  }

  lin <- select[is.na(match(select, plotvars))]
  if(length(lin)==0){

  } else if (length(lin)==1){
    cat(lin, "does not have a nonlinear effect and will not have a plot.\n")
  } else {
    cat(paste(lin,collapse=", "), "do not have nonlinear effects and will not have plots.\n")
  }

  # Initialize list of plots and create plot titles
  nz <- length(plotvars)
  nonlineffs <- vector("list", nz)
  if(length(renames)==nz){
    covnames <- renames
  } else {
    covnames <- plotvars
    if(!is.null(renames)) cat("Length of renames does not match number of nonlinear effects: ignoring renames.\n")
  }

  var.index <- match(fit$nonlinear, plotvars)
  var.index <- var.index[!is.na(var.index)]
  ylab <- names(attr(fit$terms, "dataClasses"))[1]

  nlin <- length(fit$linear)
  data <- cbind( fit$x[,1+1:nlin], fit$z)
  data[] <- 0
  N <- nrow(data)

  for(i in 1:nz){
    xlab <- covnames[i]

    data[ , nlin+var.index[i] ] <- fit$z[,var.index[i]]
    ghat <- suppressWarnings(predict(fit, newdata=data))
    mghat <- mean(ghat)
    ghat <- ghat - mghat
    z <- cbind( z=fit$z[,var.index[i]], effect=ghat )

    pp <- seq( min(fit$z[,var.index[i]]), max(fit$z[,var.index[i]]), length.out=N )
    data[ , nlin+var.index[i] ] <- pp
    ghat <- suppressWarnings(predict(fit, newdata=data))
    ghat <- ghat - mghat
    pp <- cbind(pp=pp, effect=ghat)

    nonlineffs[[i]] <- list(z=z, pp=pp, xlab=xlab, ylab=ylab, tau=fit$tau, var=plotvars[i],
                            name=xlab, call=fit$call, formula=fit$formula)

    data[ , nlin+var.index[i] ] <- 0
  }

  names(nonlineffs) <- plotvars
  class(nonlineffs) <- "plaqreffect"
  return(nonlineffs)
}



plot.plaqreffect <- function(x, select=NULL, rug = TRUE, jit = TRUE, titles = NULL,
                             pages = 0, type="l", ...)
{
  ### Thank you to Roger Koenker (quantreg) for the skeleton code.
  SetLayout <- function(m, p) {
    if (p > m) 
      p <- m
      if (p < 0) 
        p <- 0
        if (p != 0) {
          q <- m%/%p
          if ((m%%p) != 0) {
              q <- q + 1
              while (q * (p - 1) >= m) p <- p - 1
          }
          c <- trunc(sqrt(q))
          if (c < 1) 
              c <- 1
          r <- q%/%c
          if (r < 1) 
              r <- 1
          while (r * c < q) r <- r + 1
          while (r * c - q > c && r > 1) r <- r - 1
          while (r * c - q > r && c > 1) c <- c - 1
          oldpar <- par(mfrow = c(r, c))
      }
      else oldpar <- par()
      return(oldpar)
  }

  if(is.na(x)[1]){
    cat("There are no nonlinear effects to plot.\n")
    return(NULL)
  }

  y <- x
  m <- length(y)

  if (length(select)) {
    if (all(select %in% 1:m)){
      oldpar <- SetLayout(length(select), pages)
      y <- x[select]
      m <- length(y)
    }
    else stop(paste("select must be in 1:", m, sep = ""))
  }
  else oldpar <- SetLayout(m, pages)
  if ((pages == 0 && prod(par("mfrow")) < m && dev.interactive()) || 
      pages > 1 && dev.interactive()) 
      ask <- TRUE
  else ask <- FALSE
  if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }

  if(length(titles) != m){
    plot_title <- vector("list", m)
    for(j in 1:m) plot_title[[j]] <- bquote("Effect of"~.(y[[j]]$name)~(tau==.(y[[j]]$tau)))
    if(!is.null(titles)){
      cat("Length of titles does not match number of nonlinear effects: ignoring titles.\n")
    }
  } else {
    plot_title <- titles
  }

  for(j in 1:m){
    plot(y[[j]]$pp, main=plot_title[j], xlab=y[[j]]$xlab, ylab=y[[j]]$ylab, type=type,...)
    if (rug) {
      if (jit) 
        rug(jitter(y[[j]]$z[,1]))
      else rug(y[[j]]$z[,1])
    }
  }
  if (pages > 0){
    par(oldpar)
  }
}




print.plaqreffect <- function (x, ...) 
{
    if (is.na(x[1])){
      cat("There no nonlinear effects.\n")
    } else {
      m <- length(x)
      nms <- rep("NA", m)
      renms <- nms
      for(i in 1:m){
        nms[i] <- x[[i]]$var
        renms[i] <- x[[i]]$name
      } 
      cat("Call:\n")
      print(x[[1]]$call)

      if(!all(renms==nms)){
        mat <- matrix(1:m, nrow=1)
        rownames(mat) <- "Order:"
        colnames(mat) <- nms
        cat("\nNonlinear Effects:\n")
        print(mat)
        rownames(mat) <- "Order:"
        colnames(mat) <- renms
        cat("\nNonlinear Effects' Names:\n")
        print(mat)        
      } else {
          cat("\nNonlinear Effects:\n",nms,"\n")
      }
    }
}