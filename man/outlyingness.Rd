\name{outlyingness}
\alias{outlyingness}
\title{Stahel-Donoho outlyingness of points relative to a dataset}
\description{
Computes the Stahel-Donoho outlyingness (SDO) of \eqn{p}-dimensional points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}. For each multivariate point \eqn{z_i}, its outlyingness relative to \code{x} is defined as its maximal univariate Stahel-Donoho outlyingness measured over all directions. To obtain the univariate Stahel-Donoho outlyingness in the direction \eqn{v}, the dataset \code{x} is projected on \eqn{v}, and the robustly standardized distance of \eqn{v'z_i} to the robust center of the projected data points \code{x}\eqn{v} is computed.
}
\usage{
outlyingness(x, z = NULL, options = list())
}
\arguments{
  \item{x}{An \eqn{n} by \eqn{p} data matrix.}
  \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
          rowwise the points \eqn{z_i} for which to compute the outlyingness.
          If \code{z} is not specified, it is set equal to \code{x}.}
  \item{options}{ A list of available options:
                  \itemize{
                      \item{ \code{type} \cr
                      Determines the desired type of invariance and should be one
                      of \code{"Affine"}, \code{"Rotation"} or \code{"Shift"}. When
                      the option \code{"Affine"} is used,  the directions \eqn{v} are
                      orthogonal to hyperplanes spanned by \eqn{p} observations from
                      \code{x}. When the option \code{"Rotation"} is used, the directions
                      pass by two randomly selected observations from \code{x}.
                      With the option \code{"Shift"}, directions are randomly generated. \cr
                      Defaults to \code{"Affine"}.}
                      \item { \code{ndir} \cr
                      Determines the number of directions \eqn{v} by setting \code{ndir}
                      to a specific number or to \code{"all"}. In the latter case, an
                      exhaustive search over all possible directions (according to \code{type})
                      is performed. If \code{ndir} is larger than the number of possible
                      directions, the algorithm will automatically use this setting. \cr
                      Defaults to \eqn{250p} when \code{type="Affine"}, to 5000 when
                      \code{type="Rotation"} and to 12500 when \code{type="Shift"}.}
                      \item{ \code{stand} \cr
                      Determines how to robustly standardize the projected data:
                      \code{"MedMad"} uses the median and the MAD, \code{"unimcd"}
                      uses the univariate MCD of location and scale. \cr
                      Defaults to \code{"MedMad"}.}
                      \item{ \code{centered} \cr
                      When the data matrix \code{x} is already centered, no robust center
                      should be computed in each direction. In that case, \code{centered} should
                      be set to \code{TRUE}. \cr
                      Defaults to \code{FALSE}.}
                      \item{ \code{h} \cr
                      When the input argument \code{stand} is equal to \code{unimcd},
                      the parameter \code{h} controls the number of data points that
                      define the MCD estimator (see \code{\link[robustbase]{covMcd}}
                      in the \code{robustbase} package).
                      This value should lie between \eqn{[n/2]+1} and \eqn{n}. \cr
                      Defaults to \eqn{[n/2]+1}.}
                      \item{ \code{seed} } \cr
                      A strictly positive integer specifying the seed to be used
                      by the \code{C++} code.\cr
                      Defaults to \eqn{10}.
                  }
                }
}
\details{
The Stahel-Donoho outlyingness has been introduced by Stahel (1981) and Donoho (1982).
It is mostly suited to measure the degree of outlyingness of multivariate points with respect to a data cloud from an elliptical distribution.

Depending on the dimension \eqn{p}, different approximate algorithms are implemented. The affine invariant algorithm can only be used when \eqn{n > p}. It draws \code{ndir} times at random \eqn{p} observations from \code{x} and considers the direction orthogonal to the hyperplane spanned by these \eqn{p} observations. At most \eqn{p} out of \eqn{n} directions can be considered. The orthogonal invariant version can be applied to high-dimensional data. It draws \code{ndir} times at random \eqn{2} observations from \code{x} and considers the direction through these observations. Here, at most 2 out of \eqn{n} directions can be considered. Finally, the shift invariant version randomly draws \code{ndir} vectors from the unit sphere.

The resulting  Stahel-Donoho outlyingness values are invariant to affine transformations, rotations and shifts respectively provided that the \code{seed} is kept fixed at different runs of the algorithm. Note that the SDO values are guaranteed to increase when more directions are considered provided the seed is kept fixed, as this ensures that the random directions are generated in a fixed order.

An observation from \code{x} and \code{z} is flagged as an outlier if its SDO exceeds a cutoff value. This cutoff is equal to the squareroot of the 0.99 quantile of the chi-squared distribution with \eqn{p} degrees of freedom, multiplied by the median of the SDO values from the \code{x} matrix.

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it. Moreover, from the definition of the Stahel-Donoho outlyingness it follows that the outlyingness is ill-defined when the robust scale of the data projected on the direction \eqn{v} equals zero. In this case the algorithm will stop and give a warning. The returned values then include the direction \eqn{v} as well as an indicator specifying which of the observations of \code{x} belong to the hyperplane orthogonal to \eqn{v}.

}
\value{
A list with components: \cr
  \item{outlyingnessX}{Vector of length \eqn{n} giving the outlyingness of the observations in \code{x}.}
  \item{outlyingnessZ}{Vector of length \eqn{m} giving the outlyingness of the points in \code{z}.}
  \item{cutoff}{Points whose outlyingness exceeds this cutoff can be considered as outliers with respect to \code{x}.}
  \item{flagX}{Observations of \code{x} whose outlyingness exceeds the cutoff value receive a flag \code{FALSE}, regular observations receive a flag \code{TRUE}.}
  \item{flagZ}{Points of \code{z} whose outlyingness exceeds the cutoff value receive a flag equal to \code{FALSE}, otherwise they receive a flag \code{TRUE}.}
          \item{singularSubsets}{When the input parameter type is equal to \code{"Affine"},
                             the number of \eqn{p}-subsets that span a subspace of
                             dimension smaller than \eqn{p-1}. In such a case the orthogonal
                             direction can not be uniquely determined. This is an indication
                             that the data are not in general position.
                             When the input parameter type is equal to \code{"Rotation"} it
                             is possible that two randomly selected points of the data coincide
                             due to ties in the data. In such a case this value signals how many
                             times this happens.
                             }
  \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace. When a direction \eqn{v} is found such that the robust scale of \eqn{xv} is equal to zero, this equals \eqn{v}.}
  \item{inSubspace}{When a direction \eqn{v} is found such that the robust scale of \eqn{xv} is zero, the observations from \code{x} which belong to the hyperplane orthogonal to \eqn{v} receive a value TRUE. The other observations receive a value FALSE.}
}
\seealso{
  \code{\link{projdepth}}, \code{\link{projmedian}}, \code{\link{adjOutl}}, \code{\link{sprojdepth}}, \code{\link{sprojmedian}}
  }

\references{
Stahel W.A. (1981). Robuste Schatzungen: infinitesimale Optimalitat und Schatzungen  von Kovarianzmatrizen. PhD Thesis, ETH Zurich.

Donoho D.L. (1982). Breakdown properties of multivariate location estimators. Ph.D. Qualifying paper, Dept. Statistics, Harvard University, Boston.

Maronna R.A., Yohai V. (1995). The behavior of the Stahel-Donoho robust multivariate estimator. \emph{Journal of the American Statistical Association}, \bold{90}, 330--341.

}


\author{P. Segaert using \code{C++} code by K. Vakili and P. Segaert.}

\examples{

# Compute the outlyingness of a simple two-dimensional dataset.
# Outliers are plotted in red.

if (requireNamespace("robustbase", quietly = TRUE)) {
    BivData <- log(robustbase::Animals2)
} else {
  BivData <- matrix(rnorm(120), ncol = 2)
  BivData <- rbind(BivData, matrix(c(6,6, 6, -2), ncol = 2))
}

Result <- outlyingness(x = BivData)
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

# The number of directions may be specified through
# the option list. The resulting outlyingness is
# monotone increasing in the number of directions.
Result1 <- outlyingness(x = BivData,
                        options = list(ndir = 50)
                        )
Result2 <- outlyingness(x = BivData,
                        options = list(ndir = 100)
                        )
which(Result2$outlyingnessX - Result1$outlyingnessX < 0)
# This is however not the case when the seed is changed
Result1 <- outlyingness(x = BivData,
                        options = list(ndir = 50)
                        )
Result2 <- outlyingness(x = BivData,
                        options = list(ndir = 100,
                                       seed = 950)
                        )
plot(Result2$outlyingnessX - Result1$outlyingnessX,
     xlab = "Index", ylab = "Difference in outlyingness")

# We can also consider directions through two data
# points. If the sample is small enough one may opt
# to search over all choose(n,2) directions.
# Note that the computational load increases dramatically
# as n becomes larger.
Result <- outlyingness(x = BivData,
                       options = list(type = "Rotation",
                                      ndir = "all")
                       )
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

# Alternatively one may consider randomly generated directions.
Result <- outlyingness(x = BivData,
                       options = list(type = "Shift",
                                      ndir = 1000)
                      )
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

# The default option of using the MAD for the scale may be
# changed to using the univariate mcd.
Result <- outlyingness(x = BivData,
                       options = list(type = "Affine",
                                      ndir = 1000,
                                      stand = "MedMad",
                                      h = nrow(BivData))
                      )
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

}

\keyword{multivariate}
