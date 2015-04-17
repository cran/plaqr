\name{predictInt}
\alias{predictInt}
\title{Prediction Inteval for Quantile Regression} 
\description{Predicts future values using the median and finds a prediction interval for future values using an upper and lower quantile.  The lower quantile is (1-level)/2 and the upper quantile is .5 + level/2.} 
\usage{
predictInt(fit, level=.95, newdata=NULL, ...)
}
\arguments{
  \item{fit}{a fitted model of class \code{"plaqr"} or \code{"rq"} to be used for prediction.}  
  \item{level}{the prediction level required.  The lower quantile is (1-level)/2 and the upper quantile is .5 + level/2.}
  \item{newdata}{an optional data frame in which to look for variables with which 
	to predict. If omitted, the fitted values are used.}
  \item{...}{additional argument(s) for methods.}
}
\value{a matrix with columns giving the predicted median and lower and upper prediction bounds.}


\author{Adam Maidman}

\examples{
data(simData)
fit <- plaqr(y~.,~z1+z2,data=simData)
predictInt(fit, level=.95)
}

