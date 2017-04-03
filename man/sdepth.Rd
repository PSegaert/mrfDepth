\name{sdepth}
\alias{sdepth}
\title{Simplicial depth of points relative to a dataset}
\description{
Computes the simplicial depth of \eqn{p}-dimensional points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}. Only dimension \eqn{p<=2} is supported.
}
\usage{
sdepth(x, z = NULL)
}
\arguments{
	\item{x}{An \eqn{n} by \eqn{p} data matrix with observations in the rows and variables in the columns.}
  \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
           rowwise the points \eqn{z_i} for which to compute
           the simplicial depth. If \code{z} is not specified,
           it is set equal to \code{x}.}
          }
\details{
The simplicial depth has been introduced by Liu (1990). The simplicial depth of a point \eqn{z_i} is defined as the number of simplices with vertices in \code{x} that contain \eqn{z_i}. 
Exact computation of the simplicial depth for bivariate data is performed by means of the algorithm described in Rousseeuw and Ruts (1996).
To visualize the depth of bivariate data one can apply the \code{mrainbowplot} function. It plots the data with coloring according to their depth. 

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it.
}
\value{
A list with components:

  \item{depthZ}{Vector of length \eqn{m} giving the simplicial depth of the points in \code{z}. }
 \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace.}
  }
\seealso{
 \code{\link{mrainbowplot}}
}
\references{
Liu, R. (1990). On a notion of data depth based on random simplices. 
\emph{The Annals of Statistics,} \bold{18}, 405--414.

Rousseuw, P.J., Ruts, I. (1996). AS 307: Bivariate location depth. \emph{Applied Statistics, Journal of the Royal Statistical Society: Series C },  \bold{45}, 516--526.
}

\author{P. Segaert, based on Fortran code by  P.J. Rousseeuw, I. Ruts and A. Struyf.}

\examples{
data(bloodfat)
Result <- sdepth(x = bloodfat)
mrainbowplot(bloodfat, depth = Result$depthZ)
}
\keyword{multivariate}
