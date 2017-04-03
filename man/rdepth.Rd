\name{rdepth}
\alias{rdepth}
\title{Regression depth of hyperplanes}
\description{
Computes the regression depth of several hyperplanes with
 respect to a multiple regression dataset with \eqn{p} explanatory variables.
 The calculation is exact for \eqn{p <= 3}.
 An approximate algorithm is used for \eqn{p>3}.
}
\usage{
rdepth(x, z = NULL, ndir = NULL)
}
\arguments{
  \item{x}{An \eqn{n} by \eqn{p+1} regression dataset. \cr
          The first \eqn{p} columns are the explanatory variables. The last
          column corresponds to the response variable.
  }
  \item{z}{An \eqn{m} by \eqn{p+1} matrix containing rowwise the hyperplanes for which to compute the regression depth. The first column should contain the intercepts. \cr
           If \code{z} is not specified, it is set equal to \code{x}.}
	\item{ndir}{Controls the number of directions when the approximate
              algorithm is used. Defaults to \eqn{250p}.}
}
\details{
Regression depth has been introduced in Rousseeuw and Hubert (1999). To compute the regression depth of a hyperplane, different algorithms are used. When \eqn{p <= 3} it can be computed exactly. In higher dimensions an approximate algorithm is used (Rousseeuw and Struyf 1998). 

It is first checked whether the data \code{x} lie in a subspace of
dimension smaller than \eqn{p + 1}. If so, a warning is given, as well as the  dimension of the subspace and a direction which is orthogonal to it.}
 
\value{
 A  list with components: \cr
  \item{depthZ}{A vector containing the regression depth of the hyperplanes in \code{z}.}
  \item{singularSubsets}{ A vector of length \eqn{m}. \cr
               For each hyperplane, it contains the number of singular subsamples when the approximate algorithm is used. The actual number of directions used in the algorithm therefore equals \code{ndir-singularSubsets}.}
  \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace.}
}
\seealso{
  \code{\link{rdepthmedian}}, \code{\link{cmltest}}
}
\references{
Rousseeuw, P.J., Hubert, M. (1999). Regression depth.  \emph{Journal of the American Statistical Association}, \bold{94}, 388--402.

Rousseeuw, P.J., Struyf A. (1998). Computing location depth and regression depth in higher dimensions.  \emph{Statistics and Computing}, \bold{45}, 193--203.
}
\author{P. Segaert using Fortran code by P.J. Rousseeuw, I. Ruts and A. Struyf}
\examples{
# Illustrate the concept of regression depth in the case of simple 
# linear regression. 

data("stars")

# Compute the least squares fit. Due to outliers
# this fit will be bad and thus should result in a small
# regression depth.
temp <- lsfit(x = stars[,1], y = stars[,2])$coefficients
intercept <- temp[1]
slope <- temp[2]
plot(stars, ylim = c(4,7))
abline(a = intercept, b = slope)
abline(a = -9.2, b = 3.2, col = "red")

# Let us compare the regression depth of the least squares fit
# with the depth of the better fitting red regression line.
z <- rbind(cbind(intercept, slope),
           cbind(-9.2, 3.2))
result <- rdepth(x = stars, z = z)
result$depthZ
text(x = 3.8, y = 5.3, labels = round(result$depthZ[1], digits =2))
text(x = 4.45, y = 4.8, labels = round(result$depthZ[2], digits =2), col = "red")

# Compute the depth of some other regression lines to illustrate the concept.
# Note that the stars data set has 47 observations and 1/47 = 0.0212766.
z <- rbind(cbind(6.2, 0),
           cbind(6.5, 0))
result <- rdepth(x = stars, z = z)
abline(a = 6.2, b = 0, col = "blue")
abline(a = 6.5, b = 0, col = "darkgreen")
text(x = 3.8, y = 6.25, labels = round(result$depthZ[1], digits =2), col = "blue")
text(x = 4, y = 6.55, labels = round(result$depthZ[2], digits =2), col = "darkgreen")
# One point needs to removed to make the blue line a nonfit. The green line is clearly
# a nonfit. 
}
\keyword{regression}
