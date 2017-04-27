\name{projmedian}
\alias{projmedian}
\title{Location estimates based on projection depth.}
\description{
Computes a projection depth base location estimate of a \eqn{p}-dimensional 
dataset \code{x}.
}
\usage{
projmedian(x, projection.depths = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
  \item{projection.depths}{Vector containing the projection depths of the points in \code{x}.}
  \item{options}{ A list of options to pass to the \code{projdepth} routine.
                  See \code{projdepth} for more details. }

}
\details{
 The algorithm depends on the function \code{projdepth} to calculate the
 projection depths of the dataset \code{x}. If the projection depths of the
 points \eqn{x_i} have already been calculated they can be passed as an optional
 argument to save computation time. If not, projections depths will be
 calculated and the user can pass a list with options to the
 \code{projdepth} function.

 It is first checked whether the data lie in a subspace of dimension smaller 
 than \eqn{p}. If so, a warning is given, as well as the dimension of the 
 subspace and a direction which is orthogonal to it.

}
\value{
  A list with components: \cr
     \item{max}{The point of \code{x} with maximal projection depth.
                If multiple points have maximum depth, their center
                of gravity is returned}
     \item{gravity}{The center of gravity for the 50 percent points with
                    highest projection depth. Also called the
                    center of gravity of the bag.}
     \item{Huber}{The weighted center of gravity of all points.
                  The weights are defined by the Huber
                  function with parameter
                        \code{alpha = 1 / (1+sqrt(qchisq(0.95, p)))}
                  Observations for which the projection depth
                  is more than alpha receive weight 1,
                  other points receive weight
                  \code{( sqrt(qchisq(0.95, p)) / outlyingness )^2.}}
}
\seealso{
  \code{\link{outlyingness}}, \code{\link{projdepth}}, \code{\link{adjOutl}}, \code{\link{dirOutl}}
}
\references{
 Zuo Y. (2003). Projection based depth functions and associated medians. \emph{The Annals of Statistics}, \bold{31}, 1460--1490.
}
\author{P. Segaert}
\examples{

# Compute a location estimate of a simple two-dimensional dataset.
if (requireNamespace("robustbase", quietly = TRUE)) {
  BivData <- log(robustbase::Animals2)
} else {
  BivData <- matrix(rnorm(120), ncol = 2)
  BivData <- rbind(BivData, matrix(c(6,6, 6, -2), ncol = 2))
}

result <- projmedian(x = BivData)
plot(BivData)
points(result$max, col = "red", pch = 15)
points(result$gravity, col = "blue", pch = 16)
points(result$Huber, col = "orange", pch = 17)

# Options for the underlying projdepth routine may be passed 
# using the options argument. 
result <- projmedian(x = BivData,
                     options = list(type = "Affine",
                                    ndir = 10,
                                    stand = "MedMad",
                                    h = nrow(BivData)
                                    )
                    )
plot(BivData)
points(result$max, col = "red", pch = 15)
points(result$gravity, col = "blue", pch = 16)
points(result$Huber, col = "orange", pch = 17)

# One may also calculate the depths of the points in the data
# separately. This avoids having to recompute the depths when these
# are previously calculated. 
depth.result <- projdepth(x = BivData)
result <- projmedian(x = BivData, 
                     projection.depths = depth.result$depthX)

}
\keyword{multivariate}
