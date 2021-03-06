% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.PCLasso.R
\name{plot.PCLasso}
\alias{plot.PCLasso}
\title{Plot coefficients from a PCLasso object}
\usage{
\method{plot}{PCLasso}(x, norm = TRUE, ...)
}
\arguments{
\item{x}{Fitted \code{PCLasso} model.}

\item{norm}{If TRUE, plot the norm of each group, rather than the individual
coefficients.}

\item{...}{Other graphical parameters to \code{plot}.}
}
\description{
Produces a plot of the coefficient paths for a fitted
\code{PCLasso} object.
}
\examples{
# load data
data(GBM)
data(PCGroup)

x = GBM$GBM.train$Exp
y = GBM$GBM.train$survData

# fit the PCLasso model
fit1 <- PCLasso(x, y, group = PCGroup, penalty = "grLasso")

# plot the norm of each group
plot(fit1, norm = TRUE)

# plot the individual coefficients
plot(fit1, norm = FALSE)
}
\seealso{
\code{\link{PCLasso}}
}
