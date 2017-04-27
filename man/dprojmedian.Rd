\name{dprojmedian}
\alias{dprojmedian}
\title{Location estimates based on directional projection depth.}
\description{
Computes a directional projection depth based location estimate of a 
\eqn{p}-dimensional dataset \code{x}.
}
\usage{
dprojmedian(x, dprojection.depths = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{dprojection.depths}{Vector containing the directional projection
                  depths of the points in \code{x}.}
  \item{options}{ A list of options to pass to the \code{dprojdepth} routine.
                  See \code{dprojdepth} for more details. }

}
\details{
 The algorithm depends on the function \code{dprojdepth} to calculate the
 directional projection depths of the dataset \code{x}. If the 
 directional projection depths of the points \eqn{x_i} have already been
 calculated they can be passed as an optional argument to save computation time.
 If not, directional projection depths will be calculated and the user 
 can pass a list with options to the \code{dprojdepth} function.

 It is first checked whether the data lie in a subspace of dimension smaller 
 than \eqn{p}. If so, a warning is given, as well as the dimension of the 
 subspace and a direction which is orthogonal to it.
}
\value{
  A list with components: \cr
    \item{max}{The point of \code{x} with maximal directional projection 
               depth. If multiple points have maximum depth, their center
               of gravity is returned}
    \item{gravity}{The center of gravity for the 50 percent points with
                   highest directional projection depth. Also called
                   the center of gravity of the bag.}
}
\seealso{
  \code{\link{dirOutl}}, \code{\link{dprojdepth}}, \code{\link{adjOutl}}, \code{\link{outlyingness}}
}
\references{
Rousseeuw, P.J., Raymaekers, J., Hubert, M., (2016), A Measure of Directional Outlyingness with Applications to Image Data and Video. \emph{arXiv: 1608.05012}.\\
}
\author{J. Raymaekers}
\examples{
# Compute a location estimate of a simple two-dimensional dataset.
data(bloodfat)

result <- dprojmedian(x = bloodfat)
plot(bloodfat)
points(result$max, col = "red", pch = 15)
points(result$gravity, col = "blue", pch = 16)


# Options for the underlying sprojdepth routine may be passed 
# using the options argument. 
result <- dprojmedian(x = bloodfat,options = list(type = "Affine",ndir = 10))
plot(bloodfat)
points(result$max, col = "red", pch = 15)
points(result$gravity, col = "blue", pch = 16)

# One may also calculate the depths of the points in the data
# separately. This avoids having to recompute the depths when these
# are previously calculated. 
depth.result <- dprojdepth(x = bloodfat)
result <- dprojmedian(x = bloodfat,dprojection.depths = depth.result$depthX)
}
\keyword{multivariate}
