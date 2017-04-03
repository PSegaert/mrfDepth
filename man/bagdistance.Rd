\name{bagdistance}
\alias{bagdistance}
\title{Bagdistance of points relative to a dataset}
\description{
Computes the bagdistance of \eqn{p}-dimensional points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}. To compute the bagdistance of a point \eqn{z_i} the bag of \code{x} which is defined as the depth region containing the 50\% observations with largest depth. Next, the ray from the halfspace median \eqn{\theta} through \eqn{z_i} is considered and \eqn{c_z} is defined as the intersection of this ray and the boundary of the bag. The bagdistance of \eqn{z_i} to \code{x} is then given by the ratio between
the Euclidean distance of \eqn{z_i} to the halfspace median and the Euclidean distance of \eqn{c_z} to the halfspace median.
}
\usage{
bagdistance(x, z = NULL, options = list())
}
\arguments{
	\item{x}{An \eqn{n} by \eqn{p} data matrix.}
  \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
             rowwise the points \eqn{z_i} for which to compute
             the adjusted outlyingness. If \code{z} is not specified,
             it is set equal to \code{x}. is not specified,
           it is set equal to \code{x}.}
  \item{options}{A list of available options:
                \itemize{
                      \item{ \code{approx} \cr
                           In two dimension one may choose to use an approximate algorithm or the exact algorithm to find the bag. \cr
                           Defaults to \code{TRUE}.
                           }
                      \item{ \code{max.iter} \cr
                            The maximum number of steps in the bisection algorithm to find the intersection point \eqn{c_z} (see Details).\cr
                            Defaults to \eqn{100}.
                            }
                \item{}
                All options may be specified that are passed to the \code{hdepth} function, see \code{hdepth} for details. Note that the option parameter \code{approx} is by default set to \code{TRUE} to save computation time.}}
}
\value{
  A list with components: \cr
  \item{bagdistance}{The bagdistance of the points of \code{z} with respect to the data matrix \code{x}.}
  \item{cutoff}{Points of \code{z} whose bagdistance exceeds this cutoff can be considered as outliers with respect to \code{x}. }
  \item{flag}{Points of \code{z} whose bagdistance exceeds the cutoff receive a flag equal to \code{FALSE}, otherwise they receive a flag \code{TRUE}.}
  \item{converged}{Vector of length \code{m} indicating for each point of \code{z} whether the bisection algorithm converged within the maximum number of steps specified by \code{max.iter} in the \code{options} list.}
  \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace.}
}
\details{
The bagdistance has been introduced in Hubert et al. (2015). It does not assume symmetry and is affine invariant. Note that when the halfspace is not computed in an affine invariant wat, the bagdistance cannot be affine invariant either. 

The function first computes the halfspace depth and the halfspace median of \code{x}. Additional options may be passed to the \code{hdepth} routine by specifying them in the \code{option} list argument.

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it.

Depending on the dimensions different algorithms are used. For \eqn{p=1} the bagdistance is computed exactly. For \eqn{p=2} the default setting (\code{options$approx=TRUE}) uses an approximated algorithm. Exact computation, based on the exact algoritm to compute the contours of the bag (see the \code{depthContour} function), is obtained by setting \code{options$approx} to FALSE. Note that this may lead to an increase in computation time. 

For the approximated algorithm, the intersection point \eqn{c_z} is approximated by searching on each ray the point whose depth is equal to the median of the depth values of \code{x}. As the halfspace depth is monotone decreasing along the ray, a bisection algorithm is used. Starting limits are obtained by projecting the data on the direction and considering the data point with univariate depth corresponding to the median of the halfspace depths of \code{x}. By definition the multivariate depth of this point has to be lower or equal than its univariate depth. A second limit is obtained by considering the deepest location estimate. The maximum number of iterations bisecting the current search interval can be specified through the options argument \code{max.iter}. 

An observation from \code{z} is flagged as an outlier if its bagdistance exceeds a cutoff value. This cutoff is equal to the squareroot of the 0.99 quantile of the chi-squared distribution with \eqn{p} degrees of freedom.

}\seealso{
  \code{\link{depthContour}}, \code{\link{hdepth}}, \code{\link{bagplot}}
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection. \emph{Statistical Methods & Applications}, \bold{24}, 177--202.

Hubert M., Rousseeuw P.J., Segaert P. (2016). Multivariate and functional classification using depth and distance. \emph{Advances in Data Analysis and Classification}
}

\author{P. Segaert. }

\examples{
# Generate some bivariate data
nObs <-500
N <- matrix(rnorm(nObs * 2), nrow = nObs, ncol = 2)
A <- matrix(c(1,1,.5,.1), ncol = 2, nrow = 2)
X <- N%*%A

# In two dimensions we may either use the approximate
# or exact algorithm to compute the bag.
respons.exact <- bagdistance(x = X, options = list(approx = FALSE))
respons.approx <- bagdistance(x = X, options = list(approx = TRUE))
# The approximate algorithm leads to a good approximation.
plot(respons.exact$bagdistance, respons.approx$bagdistance)
abline(a = 0, b = 1)

# In Hubert et al. (2015) it was shown that for elliptical
# distributions the bagdistance^2 relates to the Mahalanobis
# distances. This may easily be illustrated.
mahDist <- mahalanobis(x = X, colMeans(X), cov(X))
plot(respons.exact$bagdistance^2, mahDist)

# Computation for the bagdistance relies on the calculation
# of halfspace depth using the hdepth function. Options for
# the hdepth routine can be passed down using the options
# arguments. Note that the affine invariance of the bagdistance
# depens on the affine invariant calculation of the halfspace
# depth. Choosing a different type for hdepth may lead to 
# unsatisfying results. 
options <-list(type = "Rotation",
               ndir = 375,
               approx = TRUE,
               seed = 78341)
respons.exact <- bagdistance(x = X, options = options)


}

\keyword{multivariate}

