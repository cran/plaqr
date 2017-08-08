\name{bic}
\alias{bic}
\title{BIC for the Partially Linear Additive Quantile Regression Model}
\description{
Returns the BIC for the partially linear additive quantile regression model from Lee, Noh, and Park (2014).
}
\usage{
bic(fit, ...)
}
\arguments{
  \item{fit}{a \code{"plaqr"} object obtained from a call to \code{plaqr}}
  \item{...}{additional parameters which will be ignored}
}
\value{
BIC value
}
\references{
Lee, E. R., Noh, H., and Park, B. U. (2014). Model selection via bayesian information
criterion for quantile regression models. Journal of the American Statistical Association
109, 216-229.
}
\examples{
data(simData)

ss <- vector("list", 2)

ss[[2]]$degree <- 3
fit1 <- plaqr(y~., nonlinVars=~z1+z2, data=simData, splinesettings=ss) 

ss[[2]]$degree <- 4
fit2 <- plaqr(y~., nonlinVars=~z1+z2, data=simData, splinesettings=ss) 

ss[[2]]$degree <- 5
fit3 <- plaqr(y~., nonlinVars=~z1+z2, data=simData, splinesettings=ss) 

bic(fit1)
bic(fit2)
bic(fit3)
}
\author{Adam Maidman}
