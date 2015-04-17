\name{plot.plaqreffect}
\alias{plot.plaqreffect}
\title{Nonlinear Effect Plot for a Partially Linear Additive Quantile Regression Model}
\description{
 Makes nonlinear effect plots for the nonlinear effects in a fit returned from the \code{nonlinEffect} function.  Note: you cannot use this function to plot a \code{"plaqr"} object.
}
\usage{
\method{plot}{plaqreffect}(x, select=NULL, rug = TRUE, jit = TRUE, titles = NULL, pages = 0, type="l", ...)
}
\arguments{
  \item{x}{a \code{plaqreffect} object returned from nonlinEffect.}
  \item{select}{vector of indices of nonlinear terms in \code{x} to be plotted, by default all.}
  \item{rug}{if TRUE, a rugplot for the x-coordinate is plotted.}
  \item{jit}{if TRUE, the x-values of the rug plot are jittered.}
  \item{titles}{title(s) as vector of character strings, by default titles are chosen for each plot as ``Effect of CovariateName (tau=tau)''.}
  \item{pages}{number of pages desired for the plots.}
  \item{type}{the type of plot that should be drawn.}
  \item{...}{additional arguments for the plotting algorithm.}
}
\examples{
data(simData)
fit <- plaqr(y~.,~z1+z2,data=simData)
eff <- nonlinEffect(fit, select=c("z1","z2"), renames=c("Length", "Height"))
eff

plot(eff)
plot(eff, select=1, col="red")
plot(eff, select=c(2,1), titles=c("Effect Z1","Effect Z2"))
plot(eff, select=1, col="red", lwd=4)

par(mfrow=c(1,2))
plot(eff)
}
\author{Adam Maidman}
