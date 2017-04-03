\name{projdepth}
\alias{projdepth}
\title{Projection depth of points relative to a dataset}
\description{
Computes the projection depth of \eqn{p}-dimensional points \code{z}
 relative to a \eqn{p}-dimensional dataset \code{x}. 
}
\usage{
projdepth(x, z = NULL, options = list())
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
    \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
             rowwise the points \eqn{z_i} for which to compute
             the projection depth. If \code{z} is not specified,
             it is set equal to \code{x}.}
  \item{options}{ A list of options to pass to the underlying \code{outlyingness} routine. \cr
                  See \code{outlyingness} for the full list of options. }
}
\details{
Projection depth is based on the Stahel-Donoho outlyingness (SDO) and is computed as \eqn{1/(1+SDO)}. It is mostly suited to measure the degree of outlyingness of multivariate points with respect to a data cloud from an elliptical distribution.

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it.

See \code{outlyingness} for more details on the computation of the SDO. To visualize the depth of bivariate data one can apply the \code{mrainbowplot} function. It plots the data colored according to their depth. 

The output values of this function are based on the output of the \code{outlyingness} function. More details can be found there. 

}
\value{
  A list with components: \cr
      \item{depthX}{Vector of length \eqn{n} giving the projection depth
                    of the observations in \code{x}.}
      \item{depthZ}{Vector of length \eqn{m} giving the projection depth
                    of the points in \code{z}.}
      \item{cutoff}{Points whose projection depth is smaller than this
                    cutoff can be considered as outliers. Equivalently the 
                    outliers are those points whose Stahel-Donoho outlyingness
                    exceed the corresponding cutoff.}
      \item{flagX}{Observations of \code{x} whose projection depth is smaller than
                   the cutoff receive a flag \code{FALSE}, regular observations receive a flag \code{TRUE}.}
      \item{flagZ}{Points of \code{z} whose projection depth is smaller than the cutoff
                   receive a flag equal to \code{FALSE}, otherwise they receive a flag \code{TRUE}.}
      \item{singularSubsets}{When the input parameter type is equal to \code{"Affine"}, the number of \eqn{p}-subsets that span a subspace of dimension smaller than \eqn{p-1}. In that case the orthogonal direction can not be uniquely determined. This is an indication that the data are not in general position.   When the input parameter type is equal to \code{"Rotation"} it is possible that two randomly selected points of the data coincide due to ties in the data. In this case this value signals how many times this happens.}
     \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
     \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace. When a direction \eqn{v} is found such that the robust scale of \eqn{xv} is equal to zero, this equals \eqn{v}.}
  \item{inSubspace}{When a direction \eqn{v} is found such that the robust scale of \eqn{xv} is zero, the observations from \code{x} which belong to the hyperplane orthogonal to \eqn{v} receive a value TRUE. The other observations receive a value FALSE.}
}
\seealso{
  \code{\link{outlyingness}}, \code{\link{projmedian}}, \code{\link{mrainbowplot}}, \code{\link{adjOutl}}, \code{\link{sprojdepth}}, \code{\link{sprojmedian}}
}
\references{
 Zuo Y. (2003). Projection-based depth functions and associated medians.  \emph{The Annals of Statistics}, \bold{31}, 1460--1490.

}
\author{P. Segaert}
\examples{

# Compute the projection depth of a simple two-dimensional dataset.
# Outliers are plotted in red.

if (requireNamespace("robustbase", quietly = TRUE)) {
    BivData <- log(robustbase::Animals2)
} else {
  BivData <- matrix(rnorm(120), ncol = 2)
  BivData <- rbind(BivData, matrix(c(6,6, 6, -2), ncol = 2))
}

Result <- projdepth(x = BivData)
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

# A multivariate rainbowplot may be obtained using mrainbowplot.
plot.options = list(legend.title = "PD")
mrainbowplot(x = BivData, 
             depths = Result$depthX, plot.options = plot.options)

# Options for the underlying outlyingness routine may be passed 
# using the options argument. 
Result <- projdepth(x = BivData, 
                    options = list(type = "Affine",
                                   ndir = 1000,
                                   stand = "MedMad",
                                   h = nrow(BivData)
                                   )
                   )

}

\keyword{multivariate}
