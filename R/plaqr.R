
plaqr <- function(formula, nonlinVars=NULL, tau=.5, data=NULL, subset,  
            weights, na.action, method = "br", model = TRUE, contrasts = NULL, 
            splinesettings=NULL, ...)
{
    if(length(tau)>1) stop("tau must be a number (not a vector) strictly between 0 and 1.")
    if(class(nonlinVars)!="formula" & !is.null(nonlinVars)){
      stop("nonlinVars must be of class \"formula\" or NULL. NULL by default.\n")
    }
    plaqrcall <- match.call()
    rqcall <- plaqrcall
    nonlinvars <- NULL
    linvars <- attr(terms(formula, data=data), "term.labels")

    if(is.null(nonlinVars)){
      int <- ifelse(attr(terms(formula, data=data),"intercept")==1, "1", "0")
      rqcall$formula <- update(formula, paste(c("~",linvars,int),
                             collapse="+"))
    } else {
      nonlinvars <- attr(terms(nonlinVars, data=data), "term.labels")
      nonlinvars <- nonlinvars[!(nonlinvars %in% all.vars(formula)[1])]
      linvars <- linvars[!(linvars %in% nonlinvars)]
      nonlinvarsbs <- rep(NA, length(nonlinvars))
      if(length(splinesettings)!=length(nonlinvarsbs) && !is.null(splinesettings))
        stop("splinesettings must either be NULL or a list with length = number of nonlinear covariates \n (see details ?plaqr)")

      bslist <- vector("list", length(nonlinvars))
      for( i in 1:length(nonlinvars) ){
        bslist[[i]] <- paste("bs(",nonlinvars[i], sep="")

        if( !is.null(splinesettings[[i]]$df) )
          bslist[[i]] <- paste(bslist[[i]], ", df=",splinesettings[[i]]$df, sep="")
        if( !is.null(splinesettings[[i]]$knots) )
          bslist[[i]] <- paste(bslist[[i]], ", knots=c(",
            toString(splinesettings[[i]]$knots), ")", sep="")
        if( !is.null(splinesettings[[i]]$degree) )
          bslist[[i]] <- paste(bslist[[i]], ", degree=",splinesettings[[i]]$degree, sep="")
        if( !is.null(splinesettings[[i]]$Boundary.knots) )
          bslist[[i]] <- paste(bslist[[i]], ", Boundary.knots=c(",
            toString(splinesettings[[i]]$Boundary.knots), ")", sep="")
    
        bslist[[i]] <- paste(bslist[[i]], ")", sep="")
      }

      nonlinvarsbs <- do.call(c, bslist)
      rqcall$formula <- update(formula, paste(c("~","1",linvars,nonlinvarsbs),
                             collapse="+"))
    }

    rqcall[[1]] <- as.name("rq")
    rqcall$nonlinVars <- rqcall$splinesettings <- NULL
    model <- eval.parent(rqcall)
    class(model) <- c("plaqr", "rq")
    model$call <- plaqrcall
    model$linear <- linvars
    model$nonlinear <- nonlinvars
    if(is.null(nonlinVars)){
      model$z <- data.frame()
    } else{
      model$z <- model.frame(nonlinVars, data=data)
    }
    if(length(model$linear)==0) model$linear <- vector("character", 0)
    if(length(model$nonlinear)==0) model$nonlinear <- vector("character", 0)

    return(model)
}



transform_plaqr <- function(formula, nonlinVars=NULL, tau=.5, data=NULL, lambda=seq(0,1,by=.05),
            confint=NULL, B=99, subset, weights, na.action, method = "br", contrasts = NULL, 
            splinesettings=NULL)
{
    if(length(tau)>1) 
      stop("tau must be a number (not a vector) strictly between 0 and 1.")
    if(class(nonlinVars)!="formula" & !is.null(nonlinVars))
      stop("nonlinVars must be of class \"formula\" or NULL. NULL by default.\n")

    nonlinvars <- NULL
    linvars <- attr(terms(formula, data=data), "term.labels")
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), 
        names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

    if(is.null(nonlinVars)){
      int <- ifelse(attr(terms(formula, data=data),"intercept")==1, "1", "0")
      mf$formula <- update(formula, paste(c("~",linvars,int),
                             collapse="+"))
    } else {
      nonlinvars <- attr(terms(nonlinVars, data=data), "term.labels")
      nonlinvars <- nonlinvars[!(nonlinvars %in% all.vars(formula)[1])]
      linvars <- linvars[!(linvars %in% nonlinvars)]
      nonlinvarsbs <- rep(NA, length(nonlinvars))
      if(length(splinesettings)!=length(nonlinvarsbs) && !is.null(splinesettings))
        stop("splinesettings must either be NULL or a list with length = number of nonlinear covariates \n (see details ?plaqr)")

      bslist <- vector("list", length(nonlinvars))
      for( i in 1:length(nonlinvars) ){
        bslist[[i]] <- paste("bs(",nonlinvars[i], sep="")

       if( !is.null(splinesettings[[i]]$df) )
          bslist[[i]] <- paste(bslist[[i]], ", df=",splinesettings[[i]]$df, sep="")
        if( !is.null(splinesettings[[i]]$knots) )
          bslist[[i]] <- paste(bslist[[i]], ", knots=c(",
            toString(splinesettings[[i]]$knots), ")", sep="")
        if( !is.null(splinesettings[[i]]$degree) )
          bslist[[i]] <- paste(bslist[[i]], ", degree=",splinesettings[[i]]$degree, sep="")
        if( !is.null(splinesettings[[i]]$Boundary.knots) )
          bslist[[i]] <- paste(bslist[[i]], ", Boundary.knots=c(",
            toString(splinesettings[[i]]$Boundary.knots), ")", sep="")

        bslist[[i]] <- paste(bslist[[i]], ")", sep="")
      }

      nonlinvarsbs <- do.call(c, bslist)
      mf$formula <- update(formula, paste(c("~","1",linvars,nonlinvarsbs),
                             collapse="+"))
    }

    mf <- eval.parent(mf)
    mt <- attr(mf, "terms")
    Y <- model.response(mf)
    X <- model.matrix(mt, mf, contrasts)

    if( min(Y) < 0 )
      stop("Response must be positive valued (add a constant)")

    ytrans <- matrix(Y, length(Y), ncol=length(lambda), byrow=FALSE)
    ytrans <- rbind(lambda, ytrans)
    ytrans <- apply( ytrans, 2, function(yy)  
      if( yy[1]==0 ){
        log(yy[-1])
      } else {
        yyy <- yy[-1]^(yy[1])
        1/(2*yy[1])*(yyy-1/yyy)
      } )

    XB <- X %*% t(  rqs.fit( X, ytrans, tau=tau )  )
    XB <- rbind(lambda, XB)

    XBbacktoY <- apply( XB, 2, function(xx)  
      if( xx[1]==0 ){
        exp(xx[-1])
      } else {
        ( xx[1]*xx[-1] + sqrt(xx[1]^2*xx[-1]^2 + 1) )^(1/xx[1])
      } )

    loss <- colSums( apply( Y - XBbacktoY, 2, function(xx) xx*(tau-1*(xx<0)) ) )

    parameter <- lambda[ which.min(loss) ]
    Y <- ytrans[ , which.min(loss) ]
    if( parameter == 1 )
      Y <- model.response(mf)

    U <- P <- confint2 <- NULL

    if( !is.null(confint) ){
      if( is.numeric(confint) ){

        if( confint >= 1 || confint <= 0 )
          confint <- .95

        U <- matrix( sample(1:length(Y),length(Y)*B, replace=TRUE), B,length(Y) )
        P <- rep(NA, B)

        for( i in 1:B ){
          Yb <- Y[ U[i,] ]
          ytransb <- ytrans[ U[i,] , ]
          Xb <- X[ U[i,] , ]
          XB <- Xb %*% t(  rqs.fit( Xb, ytransb, tau=tau )  )
          XB <- rbind(lambda, XB)

          XBbacktoY <- apply( XB, 2, function(xx)  
            if( xx[1]==0 ){
              exp(xx[-1])
            } else {
              ( xx[1]*xx[-1] + sqrt(xx[1]^2*xx[-1]^2 + 1) )^(1/xx[1])
            } )

          loss <- colSums( apply( Yb - XBbacktoY, 2, function(xx) xx*(tau-1*(xx<0)) ) )

          P[i] <- lambda[ which.min(loss) ]
        }

        confint2 <- quantile( P, c(.5-confint/2, .5+confint/2) )

      }
    }

    return( list(parameter=parameter, Y=Y, confint=confint2, U=U, P=P) )
}



trans_parameter <- function(x, parameter, inverse=FALSE){
  if( parameter == 1 ){
    x <- x


  }  else if ( inverse ){
    if( parameter == 0 ){
      x <- exp(x)
    } else {
      x <- ( parameter*x + sqrt(parameter^2*x^2 + 1) )^(1/x)
    }


  } else {
    if( min(x) < 0 )
      stop("x must be positive valued (add a constant)")
    if( parameter == 0 ){
      x <- log(x)
    } else {
      xx <- x^parameter
      x <- 1/(2*parameter)*(xx-1/xx)
      }
  }

  return(x)
}


bic <- function(fit, ...){
  n <- nrow(fit$z)
  log(  sum( fit$residuals * (fit$tau - 1*(fit$residuals<0)) )  ) + 
        length(fit$coefficients)*log(n)/(2*n)
}