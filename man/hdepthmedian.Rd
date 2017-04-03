\name{hdepthmedian}
\alias{hdepthmedian}
\title{Location estimates based on halfspace depth.}
\description{Computes the halfspace median and its corresponding
 halfspace depth for a \eqn{p}-dimensional data set \code{x}.  Computation is exact for \eqn{p \le 2} and approximate for \eqn{p > 2}.
}
\usage{
hdepthmedian(x, maxdir = NULL)
}
\arguments{
  \item{x}{An \eqn{n} by \eqn{p} data matrix.
          }
	\item{maxdir}{The number of projections used in the approximate algorithm. \cr
                Defaults to \eqn{250p}. 
                 }
}
\details{
The halfspace median, or Tukey median, is the multivariate point with largest halfspace depth with respect to the data \code{x}. This point is not always unique. In that case the halfspace median corresponds to the center of gravity of the convex set of deepest points.

It is first checked whether the data is found to lie in a subspace of
dimension lower than \eqn{p}. If so, the routine will give a warning, giving
back the dimension of the subspace together with a direction describing a
hyperplane containing this subspace.

For bivariate data the exact algorithm of Rousseeuw and Ruts (1998) is applied. 
When the data are not in general position (i.e. when there is a line containing more than two observations) dithering is performed by adding random Gaussian noise to the data. In this the ouput argument \code{dithered} will containg a flag. 

When \eqn{p > 2} the approximate algorithm of Struyf and Rousseeuw (2000) is applied. It is an iterative procedure based on projections. Their number can be chosen by the input parameter \code{maxdir}.

}
\value{
 A  list containing: \cr
  \item{median}{The coordinates of the halfspace median.\cr
            Approximate when \eqn{p>2}.}
  \item{depth}{The halfspace depth of the halfspace median. \cr
            Approximate when \eqn{p>2}.}
  \item{dithered}{Logical indicating whether dithering has been applied in the exact algorithm. \cr
                \code{FALSE} indicates no dithering has been applied. \cr
                \code{TRUE} indicates dithering has been applied.} 
	\item{ndir}{The number of projections used by the approximate algorithm. Due to the possibility of singularity of certain \eqn{p} subsamples it is possible that not all \code{maxdir} directions are evaluated. \cr
  }
	\item{AlgStopFlag}{Indicates which stopping rule is used by the approximate algorithm.\cr
                     0 indicates the maximum number of projections was reached \cr
                     1 indicates no improvement of the location estimate was
                     made after \eqn{10(p+1)} steps.
                    }
                      \item{dimension}{If the data are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{If the data are lying in a lower dimensional subspace, a direction orthogonal to this subspace.}
}
\references{
Rousseeuw P.J., Ruts I. (1998). Constructing the bivariate Tukey median. \emph{Statistica Sinica}, \bold{8}, 827--839.

Struyf A., Rousseeuw P.J. (2000). High-dimensional computation of the deepest location. \emph{Computational Statistics & Data Analysis}, \bold{34}, 415--436.
}

\author{P. Segaert based on Fortran code by P.J. Rousseeuw, I. Ruts and A. Struyf}

\examples{
# Compute a location estimate of a simple 
# two-dimensional dataset.

data(cardata90)
Result <- hdepthmedian(x=cardata90)
}
\keyword{multivariate}
