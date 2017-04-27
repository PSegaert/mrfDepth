\name{sprojmedian}
\alias{sprojmedian}
\title{Location estimates based on skewness-adjusted projection depth.}
\description{
Computes a skewness-adjusted projection depth based location estimate of a 
\eqn{p}-dimensional dataset \code{x}.
}
\usage{
sprojmedian(x, sprojection.depths = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{sprojection.depths}{Vector containing the skewness-adjusted projection
                  depths of the points in \code{x}.}
  \item{options}{ A list of options to pass to the \code{sprojdepth} routine.
                  See \code{sprojdepth} for more details. }

}
\details{
 The algorithm depends on the function \code{sprojdepth} to calculate the
 skewness-adjusted projection depths of the dataset \code{x}. If the 
 skewness-adjusted projection depths of the points \eqn{x_i} have already been
 calculated they can be passed as an optional argument to save computation time.
 If not, skewness-adjusted projection depths will be calculated and the user 
 can pass a list with options to the \code{sprojdepth} function.

 It is first checked whether the data lie in a subspace of dimension smaller 
 than \eqn{p}. If so, a warning is given, as well as the dimension of the 
 subspace and a direction which is orthogonal to it.
}
\value{
  A list with components: \cr
    \item{max}{The point of \code{x} with maximal skewness-adjusted projection 
               depth. If multiple points have maximum depth, their center
               of gravity is returned}
    \item{gravity}{The center of gravity for the 50 percent points with
                   highest skewness-adjusted projection depth. Also called
                   the center of gravity of the bag.}
}
\seealso{
  \code{\link{adjOutl}}, \code{\link{sprojdepth}}, \code{\link{dirOutl}}, \code{\link{outlyingness}}
}
\references{
 Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.
 
 Hubert M, Rousseeuw P.J., Segaert P. (2015). Multivariate Functional Outlier Detection. \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}
\author{P. Segaert}
\examples{
# Compute a location estimate of a simple two-dimensional dataset.
data(bloodfat)

result <- sprojmedian(x = bloodfat)
plot(bloodfat)
points(result$max, col = "red", pch = 15)
points(result$gravity, col = "blue", pch = 16)


# Options for the underlying sprojdepth routine may be passed 
# using the options argument. 
result <- sprojmedian(x = bloodfat,
                      options = list(type = "Affine",
                                     ndir = 10
                                     )
                     )
plot(bloodfat)
points(result$max, col = "red", pch = 15)
points(result$gravity, col = "blue", pch = 16)

# One may also calculate the depths of the points in the data
# separately. This avoids having to recompute the depths when these
# are previously calculated. 
depth.result <- sprojdepth(x = bloodfat)
result <- sprojmedian(x = bloodfat, 
                      sprojection.depths = depth.result$depthX)
}
\keyword{multivariate}
