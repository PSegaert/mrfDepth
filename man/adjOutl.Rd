\name{adjOutl}
\alias{adjOutl}
\title{Adjusted outlyingness of points relative to a dataset}
\description{
Computes the skew-adjusted outlyingness of \eqn{p}-dimensional points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}. For each multivariate point \eqn{z_i}, its adjusted outlyingness relative to \code{x} is defined as its maximal univariate adjusted outlyingness measured over all directions. To obtain the univariate adjusted outlyingness in the direction \eqn{v}, the dataset \code{x} is projected on \eqn{v}, and the robustly skew-adjusted standardized distance of \eqn{v'z_i} to the median of the projected data points \code{x}\eqn{v} is computed.
}
\usage{
adjOutl(x, z = NULL, options = list())
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
                      one of \code{"Affine"}, \code{"Rotation"} or \code{"Shift"}.
                      When the option \code{"Affine"} is used,  the directions
                      \eqn{v} are orthogonal to hyperplanes spanned by \eqn{p}
                      observations from \code{x}. When the option \code{"Rotation"}
                      is used, the directions pass by two randomly selected
                      observations from \code{x}. With the option \code{"Shift"},
                      directions are randomly generated.\cr
                      Defaults to \code{"Affine"}.}
                      \item{ \code{ndir} \cr
                      Determines the number of directions \eqn{v} by setting
                      \code{ndir} to a specific number or to \code{"all"}.
                      In the latter case, an exhaustive search over all possible
                      directions (according to \code{type}) is performed.
                      If \code{ndir} is larger than the number of possible directions,
                      the algorithm will automatically use this setting. \cr
                      Defaults to
                      \eqn{250p} when \code{type="Affine"}, to 5000 when
                      \code{type="Rotation"} and to 12500 when \code{type="Shift"}.}
                      \item{ \code{seed} } \cr
                      A strictly positive integer specifying the seed to be used
                      by the \code{C++} code.\cr
                      Defaults to \eqn{10}.
                  }
                }
}
\details{
The adjusted outlyingness (AO) of multivariate data was introduced in Brys et al. (2005) and studied in more detail in Hubert and Van der Veeken (2008). It extends the Stahel-Donoho outlyingness towards skewed distributions.

Depending on the dimension \eqn{p}, different approximate algorithms are implemented. The affine invariant algorithm can only be used when \eqn{n > p}. It draws \code{ndir} times at random \eqn{p} observations from \code{x} and considers the direction orthogonal to the hyperplane spanned by these \eqn{p} observations. At most \eqn{p} out of \eqn{n} directions can be considered. The orthogonal invariant version can be applied to high-dimensional data. It draws \code{ndir} times at random \eqn{2} observations from \code{x} and considers the direction through these two observations. Here, at most 2 out of \eqn{n} directions can be considered. Finally, the shift invariant version randomly draws \code{ndir} vectors from the unit sphere.

The resulting AO values are invariant to affine transformations, rotations and shifts respectively provided that the \code{seed} is kept fixed at different runs of the algorithm. Note that the AO values are guaranteed to increase when more directions are considered provided the seed is kept fixed, as this ensures that the random directions are generated in a fixed order.

An observation from \code{x} and \code{z} is flagged as an outlier if its AO exceeds a cutoff value. This cutoff value is equal to the squareroot of the 0.99 quantile of the chi-squared distribution with \eqn{p} degrees of freedom, multiplied by the median of the AO values from the \code{x} matrix.

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it. Furthermore, the univariate adjusted outlyingness of the projected points \eqn{xv} is ill-defined when the scale in its denominator becomes zero. This can happen when many observations collapse. In these cases the algorithm will stop and give a warning. The returned values then include the direction \eqn{v} as well as an indicator specifying which of the observations of \code{x} belong to the hyperplane orthogonal to \eqn{v}.

This function extends the \code{adjOutlyingness} function in the package \code{robustbase}. It has more options for choosing the directions, it allows to compute the adjusted outlyingness of points not belonging to the data matrix \code{x} and it is faster as it is fully implemented in C++. On the other hand, the constants (3 and -4) used in the definition of the adjusted outlyingness can not be modified in this implementation.
}
\value{
  A list with components: \cr
  \item{outlyingnessX}{Vector of length \eqn{n} giving the adjusted
    outlyingness of the observations in \code{x}.}
  \item{outlyingnessZ}{Vector of length \eqn{m} giving the adjusted
    outlyingness of the points in \code{z} relative to \code{x}.}
  \item{cutoff}{Points whose adjusted outlyingness exceeds this cutoff can be considered as outliers with respect to \code{x}. }
  \item{flagX}{Observations of \code{x} whose adjusted outlyingness exceeds the cutoff receive a flag \code{FALSE}, regular observations receive a flag \code{TRUE}.}
  \item{flagZ}{Points of \code{z} whose adjusted outlyingness exceeds the cutoff receive a flag equal to \code{FALSE}, otherwise they receive a flag \code{TRUE}.}
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
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace, a direction orthogonal to this subspace. When a direction \eqn{v} is found such that the robust skew-adjusted scale of \eqn{xv} is equal to zero, this equals \eqn{v}.}
  \item{inSubspace}{When a direction \eqn{v} is found such that AO(\eqn{xv}) is ill-defined, the observations from \code{x} which belong to the hyperplane orthogonal to \eqn{v} receive a value TRUE. The other observations receive a value FALSE.}
}
\seealso{
  \code{\link{sprojdepth}}, \code{\link{sprojmedian}}, \code{\link{dirOutl}}, \code{\link{outlyingness}}

  \code{\link[robustbase]{adjbox}}, \code{\link[robustbase]{adjOutlyingness}} from package robustbase.

}
\references{
 Brys G., Hubert M., Rousseeuw P.J. (2005). A robustification of Independent Component Analysis. \emph{Journal of Chemometrics}, \bold{19}, 364--375.

 Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.

 Hubert M., Vandervieren E. (2008). An adjusted boxplot for skewed distributions. \emph{Computational Statistics & Data Analysis}, \bold{52}, 5186--5201.

}
\author{P. Segaert using \code{C++} code by K. Vakili, P. Segaert, G. Brys and M. Maechler.}
\examples{
# Compute the adjusted outlyingness of a simple
# two-dimensional dataset. Outliers are plotted
# in red.
data("geological")
BivData <- geological[c("MnO","MgO")]
Result <- adjOutl(x = BivData)
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

# The number of directions may be specified through
# the option list. The resulting adjusted outlyingness
# is monotone increasing in the number of directions.
Result1 <- adjOutl(x = BivData,
                   options = list(ndir = 50)
                   )
Result2 <- adjOutl(x = BivData,
                   options = list(ndir = 100)
                   )
which(Result2$outlyingnessX - Result1$outlyingnessX < 0)
# This is however not the case when the seed is changed
Result1 <- adjOutl(x = BivData,
                   options = list(ndir = 50)
                   )
Result2 <- adjOutl(x = BivData,
                   options = list(ndir = 100,
                                  seed = 950)
                   )
plot(Result2$outlyingnessX - Result1$outlyingnessX,
     xlab = "Index", ylab = "Difference in AO")

# We can also consider directions through two data
# points. If the sample is small enough one may opt
# to search over all choose(n,2) directions.
# Note that the computational load increases dramatically
# as n becomes larger.
data("bloodfat")
BivData <- bloodfat[1:100,] # Consider a small toy example.
Result <- adjOutl(x = BivData,
                  options = list(type = "Rotation",
                                 ndir = "all")
                  )
IndOutliers <- which(!Result$flagX)
plot(BivData)
points(BivData[IndOutliers,], col = "red")

# Alternatively one may consider randomly generated directions.
data("bloodfat")
Result = adjOutl(x = bloodfat,
                 options = list(type = "Shift",
                                ndir = 1000)
                 )
IndOutliers <- which(!Result$flagX)
plot(bloodfat)
points(bloodfat[IndOutliers,], col = "red")
}

\keyword{multivariate}
