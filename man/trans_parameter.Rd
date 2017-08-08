
\name{trans_parameter}
\alias{trans_parameter}
\title{Transformation of the Response Variable}
\description{
Transform the response variable using the one-paremter, symmetric transformation of Geraci and Jones (2015).
}
\usage{
trans_parameter(x, parameter, inverse=FALSE)
}
\arguments{
  \item{x}{a vector of values to be transformed (the response variable)}
  \item{parameter}{a real-valued transformation parameter.  0 corresponds to the log transformation and 1 corresponds to the identity.  See Geraci and Jones (2015) for more information on the one-parameter, symmetric transformation.}
  \item{inverse}{If TRUE, the inverse transformation is done to transform the variable back to the original scale.  If FALSE, the standard transformation is computed.}
}
\value{
Returns a vector of the transformed (or back-transformed) variable.
}
\references{
Geraci, M. and Jones, M. (2015). Improved transformation-based quantile regression. Canadian Journal of Statistics 43, 118-132.

Maidman, A., Wang, L. (2017).  New Semiparametric Method for Predicting High-Cost Patients. Preprint.
}
\examples{
data(simData)
simData$Y <- exp(simData$y)

tparam <- transform_plaqr(Y~x1+x2+x3, nonlinVars=~z1+z2, data=simData) 

simData$newy <- trans_parameter(simData$Y, tparam$parameter)

fit <- plaqr(newy~x1+x2+x3, nonlinVars=~z1+z2, data=simData) 

trans_parameter( predictInt(fit), tparam$parameter, inverse=TRUE)
}
\author{Adam Maidman}
