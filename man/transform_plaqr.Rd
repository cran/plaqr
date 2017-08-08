\name{transform_plaqr}
\alias{transform_plaqr}
\title{Transformation for Partially Linear Additive Quantile Regression}
\description{
Returns the estimated transformation parameter for the one-parameter symmetric transformation (Geraci and Jones, 2015).  Confidence intervals for the transformation parameter can also be created using the bootstrap.  The response variable must be strictly positive; a constant can be added to the variable to ensure that all values are positive.
}
\usage{
transform_plaqr(formula, nonlinVars=NULL, tau=.5, data=NULL, lambda=seq(0,1,by=.05),
                confint=NULL, B=99, subset, weights, na.action, method = "br", 
                contrasts = NULL, splinesettings=NULL)
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
  \item{lambda}{
    a real-valued sequence of possible transformation parameters.  0 corresponds to the log transformation and 1 corresponds to the identity.  The transformation is symmetric so a negative transformation parameter is redundant and can be avoided.  See Geraci and Jones (2015) for more information on the one-parameter, symmetric transformation.
  }
  \item{confint}{
    a \code{confint} confidence interval for the transformation parameter will be created if \code{confint} is a number between 0 and 1 (otherwise automatically creates 95\% CI).  Otherwise, no confidence interval will be created.  The bootstrap is used to create the confidence interval.
  }
  \item{B}{
    the number of bootstrap replications for the confidence interval.  If no confidence interval is being created, this argument is ignored.
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be
    used in the fitting process.}
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
  \item{splinesettings}{
    a list of length equal to the number of nonlinear effects containing arguments to pass to the \code{bs} function for each term.  Each element of the list is either \code{NULL} or a list with named elements correpsonding to the arguments in \code{bs}.  If not \code{NULL}, the first element of splinesettings corresponds to the first nonlinear effect and so on.
  }
}
\value{
Returns the following:
  \item{parameter}{
    The transformation parameter
  }
  \item{Y}{
  The values of the transformed response
  }
  \item{confint}{
  If a confidence interval is created, this is the confidence interval for the transformation parameter.  Otherwise, \code{NULL}.  
  }
  \item{U}{
  If a confidence interval is created, a \code{B} by \code{n} matrix containing the indices used in each bootstrap sample.  Otherwise, \code{NULL}.
  }
  \item{P}{
  If a confidence interval is created, a \code{B} length vector containing the transformation parameter estimated in each bootstrap sample.  Otherwise, \code{NULL}.
  }  
}
\references{
Geraci, M. and Jones, M. (2015). Improved transformation-based quantile regression. Canadian Journal of Statistics 43, 118-132.

Maidman, A., Wang, L. (2017).  New Semiparametric Method for Predicting High-Cost Patients. Preprint.
}
\examples{
data(simData)

simData$Y <- exp(simData$y)

transform_plaqr(Y~x1+x2+x3, nonlinVars=~z1+z2, data=simData) 

transform_plaqr(Y~x1+x2+x3, nonlinVars=~z1+z2, confint=.95, data=simData) 
}



\author{Adam Maidman}