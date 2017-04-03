\name{rdepthmedian}
\alias{rdepthmedian}
\title{Hyperplane of maximal regression depth}
\description{
 Computes the deepest regression, i.e. the hyperplane with maximal
 regression depth given a regression dataset with \eqn{p} explanatory variables. The computation is exact for simple regression, and approximate otherwise.
}
\usage{
rdepthmedian(x, maxit = NULL)
}
\arguments{
	\item{x}{An \eqn{n} by \eqn{p+1} regression data matrix. \cr
          The first \eqn{p} columns are the explanatory variables. The last
          column corresponds to the response variable.}
	\item{maxit}{The maximum number of iterations. \cr
              Defaults to \eqn{100}.}
}
\details{
In simple regression the deepest regression fit can be computed exactly by considering all lines through two data points and taking the one with maximal regression depth. If several lines have the same maximal regression depth, their average is taken.

When \eqn{p > 1}, the approximate MEDSWEEP algorithm is applied (Van Aelst et al, 2002).

 It is first checked whether the data lie in a subspace of
 dimension smaller than \eqn{p + 1}. If so, a warning is given, as well as the  dimension of the subspace and a direction which is orthogonal to it.
}

\value{
 A  list with components: \cr
	\item{deepest}{A \eqn{(p+1)} vector containing the estimates of the deepest regression fit. The last \eqn{p} values are the slopes, the first value corresponds to the intercept.}
  \item{depth}{The depth of the deepest regression hyperplane.}
  \item{ndir}{ The number of performed iterations used in the medsweep algorithm. }
  \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
	\item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace.}

}
\seealso{
  \code{\link{rdepth}}, \code{\link{cmltest}}
}

\references{
Van Aelst S., Rousseeuw P.J., Hubert M., Struyf A. (2002). The deepest regression method. \emph{Journal of Multivariate Analysis}, \bold{81}, 138--166.
}

\author{P. Segaert using Fortran code by S. Van Aelst}

\examples{
# Illustrate the concept of deepest regression line in the case of simple 
# linear regression. 
data("stars")
plot(stars)
result <- rdepthmedian(x = stars)
abline(result$deepest)

x <- matrix(rnorm(3000), ncol = 3) + 10
rdepthmedian(x = x)
}
\keyword{regression}
