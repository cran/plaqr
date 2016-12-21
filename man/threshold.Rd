
\name{threshold}
\alias{threshold}
\title{
Classifying a Numerical Response Using a Threshold
}
\description{
Classification of a numerical response into a ``high'' class and ``low'' class using a threshold.  This function can be used with any model that has a numerical outcome and allows for prediction using the \code{predict} function.
}
\usage{
threshold(fit, t, newdata=NULL, ...)
}
\arguments{
  \item{fit}{any model with a numerical response.}
  \item{t}{the desired threshold value.  All values above \code{t} will be labeled ``1''and all values below \code{t} will be labeled ``0''.}
  \item{newdata}{an optional data frame in which to look for variables with which to predict. If omitted, no prediction is done.}
  \item{...}{additional argument(s) for methods in the \code{predict} function.}
}
\value{
  \item{pred.class}{if \code{newdata} is not \code{NULL}, then \code{pred.class} is a vector of predicted classes for \code{newdata}.  If \code{newdata} is \code{NULL}, then \code{pred.class} is \code{NULL}.}
  \item{t}{the threshold.}
  \item{train.class}{a vector of the predicted classes of the data used in \code{fit}.}
  \item{true.class}{a vector of the true classes of the data used in \code{fit}.}
  \item{train.error}{a scalar equal to the \code{mean(train.class != true.class)}.}
  \item{true.high}{the number of observations in class``1'' using the data used in \code{fit}.}
  \item{true.low}{the number of observations in class ``1'' using the data used in \code{fit}.}
  \item{false.high}{the number of observations truly in class ``0'', but predicted to be in class ``1'' using the data used in \code{fit}.}
  \item{false.low}{the number of observations truly in class ``1'', but predicted to be in class ``1'' using the data used in \code{fit}.}
  \item{call}{the \code{call} of \code{fit}.}
  \item{formula}{the formula used in \code{fit}.}
}
\examples{
data(simData)
fit <- plaqr(y~.,~z1+z2,data=simData)
testdata <- .5*simData[4,2:6]
trh <- threshold(fit, t=9, newdata=testdata)
trh$pred.class
trh
}
\author{Adam Maidman}
