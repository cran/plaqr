\name{plaqr}
\alias{plaqr}
\title{Partially Linear Additive Quantile Regression}
\description{
Returns an object of class \code{"plaqr"} and \code{"rq"} that represents a quantile regression fit.  A nonlinear term z is transformed using \code{bs(z, degree=basis)} before fitting the model.  The formula of the model (as it appears in \code{R}) becomes \code{y~ x1 + x2 + bs(z1, degree=basis) + bs(z2, degree=basis)} where \code{bs(z1, degree=basis)} is a B-spline.
}
\usage{
plaqr(formula, nonlinVars=NULL, tau=.5, data=NULL, subset, basis=3,  
            weights, na.action, method = "br", model = TRUE, contrasts = NULL, ...)
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a \code{~} operator, 
    and the linear terms, separated by \code{+} operators, on the right.  Any terms on the right         of the \code{~} operator that also appear in \code{nonlinVars} will be included in the model as spline terms, not linear terms.
  }
  \item{nonlinVars}{
    a one-sided formula object, with a \code{~} operator to the left of the nonlinear terms seperated by \code{+} operators.  A term appearing in both \code{formula} and \code{nonlinVars} will be treated as a nonlinear term.  If \code{nonlinVars} is not \code{NULL}, then an intercept will automatically be included in the model (despite a \code{-1} or \code{0} term included in \code{formula}).
  }
  \item{tau}{
    the quantile to be estimated, this is a number strictly between 0 and 1 (for now).
  }
  \item{data}{
    a data.frame in which to interpret the variables 
    named in the formula, or in the subset and the weights argument. 
    If this is missing, then the variables in the formula should be on the 
    search list.  This may also be a single number to handle some special  
    cases -- see below for details.   
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{basis}{
    degree of the piecewise polynomial to generate the B-spline basis matrix for the nonlinear terms.  Nonlinear terms are transformed using \code{bs(z, degree=basis)}.
  }
  \item{weights}{
    vector of observation weights; if supplied, the algorithm fits
    to minimize the sum of the weights multiplied into the
    absolute residuals. The length of weights must be the same as
    the number of observations.  The weights must be nonnegative
    and it is strongly recommended that they be strictly positive,
    since zero weights are ambiguous. 
  }
  \item{na.action}{
    a function to filter missing data. 
    This is applied to the model.frame after any subset argument has been used. 
    The default (with \code{na.fail}) is to create an error if any missing values are  
    found.  A possible alternative is \code{na.omit}, which 
    deletes observations that contain one or more missing values. 
  }
  \item{model}{if TRUE then the model frame is returned.  This is
    essential if one wants to call summary subsequently.
  }
  \item{method}{
    the algorithmic method used to compute the fit.  There are several
    options:   The default method is the modified  version of the
    Barrodale and Roberts algorithm for \eqn{l_1}{l1}-regression,
    used by \code{l1fit} in S, and is described in detail in 
    Koenker and d'Orey(1987, 1994),  default = \code{"br"}. 
    This is quite efficient for problems up to several thousand observations, 
    and may be used to compute the full quantile regression process.  It 
    also implements a scheme for computing confidence intervals for 
    the estimated parameters, based on inversion of a rank test described 
    in Koenker(1994).  For larger problems it is advantagous to use 
    the Frisch--Newton interior point method \code{"fn"}. 
    And very large problems one can use the Frisch--Newton approach after 
    preprocessing \code{"pfn"}.  Both of the latter methods are
    described in detail in Portnoy and Koenker(1997).   
    There is a fifth option \code{"fnc"} that enables the user to specify
    linear inequality constraints on the fitted coefficients; in this
    case one needs to specify the matrix \code{R} and the vector \code{r}
    representing the constraints in the form \eqn{Rb \geq r}.  See the
    examples.  Finally, there are two penalized methods:  \code{"lasso"}
    and \code{"scad"} that implement the lasso penalty and Fan and Li's
    smoothly clipped absolute deviation penalty, respectively.  These
    methods should probably be regarded as experimental.
  }
  \item{contrasts}{
    a list giving contrasts for some or all of the factors 
    default = \code{NULL} appearing in the model formula. 
    The elements of the list should have the same name as the variable 
    and should be either a contrast matrix (specifically, any full-rank 
    matrix with as many rows as there are levels in the factor), 
    or else a function to compute such a matrix given the number of levels. 
  }
  \item{...}{
    additional arguments for the fitting routines 
    (see the \code{rq} function in the `quantreg' package ). 
  }
}
\value{
Returns the following:
  \item{coefficients}{
    Coefficients from the fitted model
  }
  \item{x}{
  optionally the model matrix, if \code{x=TRUE}.
  }
  \item{y}{
  optionally the response, if \code{y=TRUE}.
  }
  \item{residuals}{
  the residuals from the fit.
  }
  \item{dual}{
  the vector dual variables from the fit.
  }
  \item{fitted.values}{
  fitted values from the fit.
  }
  \item{formula}{
  the formula that was used in the \code{rq} function.
  }
  \item{rho}{
  the value of the objective function at the solution.
  }
  \item{model}{
  optionally the model frame, if \code{model=TRUE}
  }
  \item{linear}{
  the linear terms used in the model fit.
  }
  \item{nonlinear}{
  the nonlinear terms used in the model fit.
  }
  \item{z}{
  the values of the nonlinear terms.
  }
  
}
\references{
Hastie, T. J. (1992) Generalized additive models. Chapter 7 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

Koenker, R. W. (2005). Quantile Regression, Cambridge U. Press.

Maidman, A., Wang, L. (2015).  Predicting High Spending Customers Using Semiparametric Quantile Regression. Preprint.

Sherwood, Benjamin Stanley. (2014). Quantile regression model selection. Retrieved from the University of Minnesota Digital Conservancy, http://hdl.handle.net/11299/163910.
}
\examples{
data(simData)
plaqr(y~., nonlinVars=~z1+z2, data=simData) 
#same as plaqr(formula= y~x1+x2+x3, nonlinVars=~z1+z2, data=simData)

plaqr(y~0, nonlinVars=~z1+z2, data=simData) #no linear terms in the model

plaqr(y~., data=simData) #all linear terms
}



\author{Adam Maidman}