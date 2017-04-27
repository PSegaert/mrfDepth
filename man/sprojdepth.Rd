\name{sprojdepth}
\alias{sprojdepth}
\title{Skewness-adjusted projection depth of points relative to a dataset}
\description{
 Computes the skewness-adjusted projection depth of \eqn{p}-dimensional
 points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}.
}
\usage{
sprojdepth(x, z = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
    \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
             rowwise the points \eqn{z_i} for which to compute
             the projection depth. If \code{z} is not specified,
             it is set equal to \code{x}.}
  \item{options}{ A list of options to pass to the underlying \code{adjOutl} routine.
                  See \code{adjOutl} for the full list of options. }
}
\details{
 Skewness-adjusted projection depth is based on the adjusted
 outlyingness and is computed as \eqn{1/(1+AO)}. As adjusted 
 outlyingness extends the Stahel-Donoho outlyingness towards 
 skewed distributions, the skewness-adjusted projection depth 
 is suited for both elliptical distributions and skewed 
 multivariate data.

 It is first checked whether the data is found to lie in a subspace of
 dimension lower than \eqn{p}. If so, a warning is given, as well as the
 dimension of the subspace and a direction which is orthogonal to it.
 
 See \code{adjOutl} for more details on the computation of the AO.
 To visualize the depth of bivariate data one can apply the 
 \code{mrainbowplot} function. It plots the data colored according to 
 their depth. 

 The output values of this function are based on the output of the 
 \code{adjOutl} function. More details can be found there. 
}
\value{
  A list with components: \cr
      \item{depthX}{Vector of length \eqn{n} giving the
                    skewness-adjusted projection depth of
                    the observations in x.}
      \item{depthZ}{Vector of length \eqn{m} giving the
                    skewness-adjusted projection depth of
                    the points in z.}
      \item{cutoff}{Points whose skew-adjusted projection depth is smaller
                than this cutoff can be considered as outliers.}
      \item{flagX}{Observations of \code{x} whose adjusted outlyingness exceeds 
                   the cutoff receive a flag \code{FALSE}, regular observations 
                   receive a flag \code{TRUE}.}
      \item{flagZ}{Points of \code{z} whose adjusted outlyingness exceeds the 
                   cutoff receive a flag equal to \code{FALSE}, otherwise they 
                   receive a flag \code{TRUE}.}
      \item{singularSubsets}{When the input parameter type is equal to \code{"Affine"},
                             the number of \eqn{p}-subsets that span a subspace of
                             dimension smaller than \eqn{p-1}. In that case the orthogonal
                             direction can not be uniquely determined. This is an indication
                             that the data are not in general position.
                             When the input parameter type is equal to \code{"Rotation"} it
                             is possible that two randomly selected points of the data coincide
                             due to ties in the data. In this case this value signals how many
                             times this is the case.
                             }
  \item{dimension}{When the data \code{x} are lying in a
                   lower dimensional subspace, the dimension
                   of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional 
                    subspace, a direction orthogonal to this subspace. When a 
                    direction \eqn{v} is found such that the robust skew-adjusted
                    scale of \eqn{xv} is equal to zero, this equals \eqn{v}.}
  \item{inSubspace}{When a direction \eqn{v} is found such that AO(\eqn{xv}) is 
                    ill-defined, the observations from \code{x} which belong to 
                    the hyperplane orthogonal to \eqn{v} receive a value TRUE. 
                    The other observations receive a value FALSE.}
}
\seealso{
  \code{\link{adjOutl}}, \code{\link{sprojmedian}}, \code{\link{mrainbowplot}}, \code{\link{dirOutl}}, \code{\link{outlyingness}}
}
\references{
 Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.
 
 Hubert M, Rousseeuw P.J., Segaert P. (2015). Multivariate Functional Outlier Detection. \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}
\author{P. Segaert with original code from M. Maechler, G. Brys, K. Vakili}
\examples{
# Compute the skewness-adjusted projection depth 
# of a simple two-dimensional dataset.
# Outliers are plotted in red.

data(bloodfat)
Result <- sprojdepth(x = bloodfat)
IndOutliers <- which(!Result$flagX)
plot(bloodfat)
points(bloodfat[IndOutliers,], col = "red")

# A multivariate rainbowplot may be obtained using mrainbowplot.
plot.options = list(legend.title = "SPD")
mrainbowplot(x = bloodfat, 
             depths = Result$depthX, plot.options = plot.options)

# Options for the underlying outlyingness routine may be passed 
# using the options argument. 
Result <- sprojdepth(x = bloodfat, 
                     options = list(type = "Affine",
                                    ndir = 1000,
                                    seed = 12345
                                    )
                    )
}
\keyword{multivariate}
