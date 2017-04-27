\name{dprojdepth}
\alias{dprojdepth}
\title{Directional projection depth of points relative to a dataset}
\description{
 Computes the directional projection depth of \eqn{p}-dimensional
 points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}.
}
\usage{
dprojdepth(x, z = NULL, options = NULL)
}
\arguments{
    \item{x}{An \eqn{n} by \eqn{p} data matrix with observations
             in the rows and variables in the columns.}
    \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
             rowwise the points \eqn{z_i} for which to compute
             the projection depth. If \code{z} is not specified,
             it is set equal to \code{x}.}
  \item{options}{ A list of options to pass to the underlying \code{dirOutl} routine.
                  See \code{dirOutl} for the full list of options. }
}
\details{
 Directional projection depth is based on the directional
 outlyingness and is computed as \eqn{1/(1+DO)}. As directional 
 outlyingness extends the Stahel-Donoho outlyingness towards 
 skewed distributions, the directional projection depth 
 is suited for both elliptical distributions and skewed 
 multivariate data.

 It is first checked whether the data is found to lie in a subspace of
 dimension lower than \eqn{p}. If so, a warning is given, as well as the
 dimension of the subspace and a direction which is orthogonal to it.
 
 See \code{dirOutl} for more details on the computation of the DO.
 To visualize the depth of bivariate data one can apply the 
 \code{mrainbowplot} function. It plots the data colored according to 
 their depth. 

 The output values of this function are based on the output of the 
 \code{dirOutl} function. More details can be found there. 
}
\value{
  A list with components: \cr
      \item{depthX}{Vector of length \eqn{n} giving the
                    directional projection depth of
                    the observations in x.}
      \item{depthZ}{Vector of length \eqn{m} giving the
                    directional projection depth of
                    the points in z.}
      \item{cutoff}{Points whose directional projection depth is smaller
                than this cutoff can be considered as outliers.}
      \item{flagX}{Observations of \code{x} whose directional outlyingness exceeds 
                   the cutoff receive a flag \code{FALSE}, regular observations 
                   receive a flag \code{TRUE}.}
      \item{flagZ}{Points of \code{z} whose directional outlyingness exceeds the 
                   cutoff receive a flag equal to \code{FALSE}, otherwise they 
                   receive a flag \code{TRUE}.}
      \item{singularSubsets}{When the input parameter type is equal to \code{"Affine"},
                             the number of \eqn{p}-subsets that span a subspace of
                             dimension smaller than \eqn{p-1}. In such a case the orthogonal
                             direction can not be uniquely determined. This is an indication
                             that the data are not in general position.
                             }
  \item{dimension}{When the data \code{x} are lying in a
                   lower dimensional subspace, the dimension
                   of this subspace.}
  \item{hyperplane}{When the data \code{x} are lying in a lower dimensional 
                    subspace, a direction orthogonal to this subspace. When a 
                    direction \eqn{v} is found such that the robust directional
                    scale of \eqn{xv} is equal to zero, this equals \eqn{v}.}
  \item{inSubspace}{When a direction \eqn{v} is found such that DO(\eqn{xv}) is 
                    ill-defined, the observations from \code{x} which belong to 
                    the hyperplane orthogonal to \eqn{v} receive a value TRUE. 
                    The other observations receive a value FALSE.}
}
\seealso{
  \code{\link{dirOutl}}, \code{\link{dprojmedian}}, \code{\link{mrainbowplot}}, \code{\link{adjOutl}}, \code{\link{outlyingness}}
}
\references{
Rousseeuw, P.J., Raymaekers, J., Hubert, M., (2016), A Measure of Directional Outlyingness with Applications to Image Data and Video. \emph{arXiv: 1608.05012}.\\
}
\author{J. Raymaekers}
\examples{
# Compute the directional projection depth 
# of a simple two-dimensional dataset.
# Outliers are plotted in red.

data(bloodfat)
Result <- dprojdepth(x = bloodfat)
IndOutliers <- which(!Result$flagX)
plot(bloodfat)
points(bloodfat[IndOutliers,], col = "red")

# A multivariate rainbowplot may be obtained using mrainbowplot.
plot.options = list(legend.title = "DPD")
mrainbowplot(x = bloodfat, 
             depths = Result$depthX, plot.options = plot.options)

# Options for the underlying outlyingness routine may be passed 
# using the options argument. 
Result <- dprojdepth(x = bloodfat,options = list(type = "Affine",ndir=100))

}
\keyword{multivariate}
