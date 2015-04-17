\name{nonlinEffect}
\alias{nonlinEffect}
\title{Nonlinear Effects Plots}
\description{
Returns an object of class \code{"plaqreffect"} which represents the effect plot(s) of the nonlinear term(s) of a \code{"plaqr"} object from the \code{plaqr} function.  A \code{"plaqreffect"} object should be plotted using the \code{plot} function.
}
\usage{
nonlinEffect(fit, select=NULL, renames=NULL) 
}
\arguments{
  \item{fit}{a \code{"plaqr"} object.}
  \item{select}{a character vector with entries matching nonlinear terms in \code{fit}.}
  \item{renames}{a character vector with length equal to the number of nonlinear terms in \code{select} (if \code{select} is \code{NULL}, the length must be equal to the number of nonlinear terms in \code{fit}).  The first entry renames the first nonlinear term for plotting purposes, and so on.  Note that \code{select} can reorder the nonlinear terms (see the examples).}
}
\value{
  A returned \code{"plaqreffect"} object to be used with the \code{"plot"} function.  Each nonlinear term is associated with a list containing information for plotting.  See the examples for accessing the list.
} 

\examples{
data(simData)
fit <- plaqr(y~.,~z1+z2,data=simData)

eff1 <- nonlinEffect(fit)
eff1
plot(eff1)

eff2 <- nonlinEffect(fit, select=c("z1","z2"), renames=c("Length", "Height"))
eff2
plot(eff2)

eff3 <- nonlinEffect(fit, select=c("z2","z1"), renames=c("Height", "Length"))
eff3
eff3$z1
eff3$z2
plot(eff3)

par(mfrow=c(1,2))
plot(eff3)
}
\author{Adam Maidman}
