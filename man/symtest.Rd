\encoding{UTF-8}
\name{symtest}
\alias{symtest}
\title{Test for angular symmetry around a specified center for bivariate data}
\description{A test based on halfspace depth for angular symmetry of the bivariate data set \code{x} around the point \code{z}. The test is only valid when \code{x} contains no duplicates.
}
\usage{
symtest(x, z, options=list())
}
\arguments{
  \item{x}{ An \eqn{n} by 2 data matrix.}
  \item{z}{ A vector of length \eqn{2}.}
\item{options}{A list of options to pass to \code{hdepth}. \cr
                    See \code{hdepth} for more information. }
}
\details{
 The following hypothesis test is performed: \cr
    H0: The data come from a continuous distribution \eqn{P} which is angularly symmetric about \code{z}. \cr
The test statistic being used is the halfspace depth of \code{z} with respect to \code{x}. Under the null hypothesis the halfspace depth of \code{z} equals \eqn{1/2}. The distribution of the teststatistic under H0 is \eqn{F_n(k) = P(hdepth(P, z) \le k)}. The \eqn{p}-value of the test is \eqn{F_n(k)} with \eqn{k= hdepth(x,z)}. 
}
\value{
	\item{pval}{The \eqn{p}-value of the hypothesis test.}
}
\references{
Rousseeuw, P.J, Struyf, A. (2002). A Depth Test for Symmetry, in:  \emph{ Goodness-of-Fit Tests and Model Validity}, \enc{Birkhäuser}{Birkhauser} Boston,  pp.401--412.
}

\author{P. Segaert}

\examples{
# Perform the test on a simple data example.
data(cardata90)
symtest(x = cardata90, z = hdepthmedian(cardata90)$median)
plot(cardata90)
}
\keyword{multivariate}
