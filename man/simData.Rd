\name{simData}
\alias{simData}
\title{Simulated Data}
\description{
A simulated data set to illustrate the functions in this package.

\code{set.seed(4)}

\code{x1  <- rbinom(100, 1,.5)}

\code{x2  <- rnorm(100)}

\code{x3  <- rnorm(100)}

\code{z1  <- runif(100,  0, 1)}

\code{z2  <- runif(100, -1, 1)}

\code{y   <- 3*x1 +1.5*x2 + 2*x3 + 5*sin(2*pi*z1) + 5*z2^3 + rnorm(100)}

\code{simData <- data.frame(y,x1,x2,x3,z1,z2)}
}
\usage{data(simData)}
\format{
  A data frame with 100 observations on the following 6 variables.
  \describe{
    \item{y}{response: expenditure}
    \item{x1}{male/female (a linear term)}
    \item{x2}{distance north/south from center (a linear term)}
    \item{x3}{distance east/west from center (a linear term)}
    \item{z1}{income/(max income) (a nonlinear term)}
    \item{z2}{spending habits on a -1 to 1 scale (frugal to lavish) (a nonlinear term)}
  }
}
