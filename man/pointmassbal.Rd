% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pointmassbal.R
\name{pointmassbal}
\alias{pointmassbal}
\title{Point Estimate Mass Balance}
\usage{
pointmassbal(y)
}
\arguments{
\item{y}{A list of matrices of observed mass flow rates. Each matrix is a separate sample component.  The rows of each matrix index the sampling location, and the columns index the sample set number.  Necessary to format exactly the same as the output of \code{twonodeSim()$simulation}.  Any arguments to \code{\link{twonodeSim}} can be used.}
}
\value{
Returns a list of vectors with mass flow rates \code{yhat} and grades \code{ahat}.  Similar format to argument \code{y}.  The index of a vector in the output is equivalent to the index of a row in \code{y}.
}
\description{
Function conducting a two component, two node, point estimate mass balance from \insertCite{willsbook}{BayesMassBal} on a two node process.  This function is provided for the purpose of comparing the performance of a Bayesian mass balance to a point estimate mass balance using any output of \code{\link{twonodeSim}}.
}
\examples{
y <- twonodeSim()$simulation

yhat <- pointmassbal(y)$yhat

}
\references{
\insertRef{willsbook}{BayesMassBal}
}
