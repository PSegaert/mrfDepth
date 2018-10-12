\name{dirOutl}
\alias{dirOutl}
\title{Directional outlyingness of points relative to a dataset}
\description{
Computes the directional outlyingness of \eqn{p}-dimensional points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}. For each multivariate point \eqn{z_i}, its directional outlyingness relative to \code{x} is defined as its maximal univariate directional outlyingness measured over all directions. To obtain the univariate directional outlyingness in the direction \eqn{v}, the dataset \code{x} is projected on \eqn{v}, and the robustly skew-adjusted standardized distance of \eqn{v'z_i} to the median of the projected data points \code{x}\eqn{v} is computed. This is done through the estimation of 2 scales, one on each side of the median, using a 1-step M-estimator of scale.
}
\usage{
dirOutl(x, z = NULL, options = list())
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix.}
    \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
             rowwise the points \eqn{z_i} for which to compute
             the adjusted outlyingness. If \code{z} is not specified,
             it is set equal to \code{x}.}
    \item{options}{ A list of available options:
                  \itemize{
                      \item{ \code{type} \cr
                      Determines the desired type of invariance and should be
                      one of \code{"Affine"}, \code{"compWise"}.
                      When the option \code{"Affine"} is used,  the directions
                      \eqn{v} are orthogonal to hyperplanes spanned by \eqn{p}
                      observations from \code{x}. With the option \code{"compWise"},
                      the directional outlyingness is computed in the directions of the
                      coordinate axes and combined through the Euclidean norm.\cr
                      Defaults to \code{"Affine"}.}
                      \item{ \code{ndir} \cr
                      When \code{type} is chosen to be \code{"Affine"}, determines the number of                              directions \eqn{v} by setting \code{ndir} to a specific number.
                      Defaults to \eqn{250p}.}
                      \item{ \code{seed}  \cr
                      A strictly positive integer specifying the seed to be used
                      to select the directions.\cr
                      Defaults to \eqn{10}.}
                  }
                }
}

\details{
The directional outlyingness (DO) of multivariate data was introduced in Rousseeuw et al. (2018). It extends the Stahel-Donoho outlyingness towards skewed distributions.

Depending on the dimension \eqn{p}, different approximate algorithms are implemented. The affine invariant algorithm can only be used when \eqn{n > p}. It draws \code{ndir} times at random \eqn{p} observations from \code{x} and considers the direction orthogonal to the hyperplane spanned by these \eqn{p} observations. At most \eqn{p} out of \eqn{n} directions can be considered. The orthogonal invariant version can be applied to high-dimensional data. It draws \code{ndir} times at random \eqn{2} observations from \code{x} and considers the direction through these two observations. Here, at most 2 out of \eqn{n} directions can be considered. Finally, the shift invariant version randomly draws \code{ndir} vectors from the unit sphere.

The resulting DO values are invariant to affine transformations, rotations and shifts respectively provided that the \code{seed} is kept fixed at different runs of the algorithm. Note that the DO values are guaranteed to increase when more directions are considered provided the seed is kept fixed, as this ensures that the random directions are generated in a fixed order.

An observation from \code{x} and \code{z} is flagged as an outlier if its DO exceeds a cutoff value. This cutoff value is determined using the procedure in Rousseeuw et al. (2018). First, the logarithm of the DO values is taken to render their distribution more symmetric, after which a normal approximation yields a cutoff on these values. The cutoff is then transformed back by applying the exponential function.

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it. Furthermore, the univariate adjusted outlyingness of the projected points \eqn{xv} is ill-defined when the scale in its denominator becomes zero. This can happen when many observations collapse. In these cases the algorithm will stop and give a warning. The returned values then include the direction \eqn{v} as well as an indicator specifying which of the observations of \code{x} belong to the hyperplane orthogonal to \eqn{v}.

}
\value{
  A list with components: \cr
  \item{outlyingnessX}{Vector of length \eqn{n} giving the directional
    outlyingness of the observations in \code{x}.}
  \item{outlyingnessZ}{Vector of length \eqn{m} giving the directional
    outlyingness of the points in \code{z} relative to \code{x}.}
  \item{cutoff}{Points whose directional outlyingness exceeds this cutoff can be considered as outliers with respect to \code{x}.}
  \item{flagX}{Observations of \code{x} whose directional outlyingness exceeds the cutoff receive a flag \code{FALSE}, regular observations receive a flag \code{TRUE}.}
  \item{flagZ}{Points of \code{z} whose directional outlyingness exceeds the cutoff receive a flag equal to \code{FALSE}, otherwise they receive a flag \code{TRUE}.}
        \item{singularSubsets}{When the input parameter type is equal to \code{"Affine"},
                             the number of \eqn{p}-subsets that span a subspace of
                             dimension smaller than \eqn{p-1}. In such a case the orthogonal
                             direction can not be uniquely determined. This is an indication
                             that the data are not in general position.
                             }
  \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace. When a direction \eqn{v} is found such that the robust skew-adjusted scale of \eqn{xv} is equal to zero, this equals \eqn{v}.}
  \item{inSubspace}{When a direction \eqn{v} is found such that DO(\eqn{xv}) is ill-defined, the observations from \code{x} which belong to the hyperplane orthogonal to \eqn{v} receive a value TRUE. The other observations receive a value FALSE.}
}
\seealso{
  \code{\link{dprojdepth}}, \code{\link{dprojmedian}}, \code{\link{adjOutl}}, \code{\link{outlyingness}}

}
\references{
Rousseeuw, P.J., Raymaekers, J., Hubert, M., (2018), A Measure of Directional Outlyingness with Applications to Image Data and Video. \emph{Journal of Computational and Graphical Statistics}, \bold{27}, 345--359.
}
\author{J. Raymaekers and P. Roussseeuw}
\examples{
# Compute the directional outlyingness of a simple
# two-dimensional dataset. Outliers are plotted
# in red.
data("geological")
BivData <- geological[c("MnO","MgO")]
Result <- dirOutl(x = BivData)
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

# The number of directions may be specified through
# the option list. The resulting adjusted outlyingness
# is monotone increasing in the number of directions.
Result1 <- dirOutl(x = BivData,options = list(ndir = 50))
Result2 <- dirOutl(x = BivData,options = list(ndir = 100))
which(Result2$outlyingnessX - Result1$outlyingnessX < 0)
# This is however not the case when the seed is changed
Result1 <- dirOutl(x = BivData,options = list(ndir = 50))
Result2 <- dirOutl(x = BivData,options = list(ndir = 100,seed = 950))

plot(Result2$outlyingnessX - Result1$outlyingnessX,
     xlab = "Index", ylab = "Difference in DO")

# Consider another example:

data("bloodfat")
BivData <- bloodfat[1:100,] # Consider a small toy example.
Result <- dirOutl(x = BivData,options = list(type = "Affine"))
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

}

\keyword{multivariate}
