\name{bagplot}
\alias{bagplot}
\title{Draws a bagplot, a bivariate boxplot}
\description{
This function draws a bagplot of bivariate data, based on the result of a call to \code{compBagplot}. The bagplot is a generalisation of the univariate boxplot to bivariate data. It aims to visualize the location, spread, skewness and outliers of the data set. 
}
\usage{
bagplot(compbag.result,
        colorbag = NULL, colorloop = NULL, colorchull = NULL,
        databag = TRUE, dataloop = TRUE, plot.fence = FALSE)
}
\arguments{
  \item{compbag.result}{The return of a call to \code{compBagplot}.}
  \item{colorbag}{The color of the bag (which contains the 50\% observations with largest depth).}
	\item{colorloop}{The color of the loop (which contains the regular observations).}
  \item{colorchull}{When the bagplot is based on halfspace depth, the depth region with maximal depth is plotted. This argument controls its color.}
  \item{databag}{Logical indicating whether data points inside the bag need to be plotted. \cr
  Defaults to \code{TRUE}.}
  \item{dataloop}{Logical indicating whether data points inside the fence need to be plotted. \cr
  Defaults to \code{TRUE}.}
    \item{plot.fence}{Logical indicating whether the fence should be plotted. \cr
  Defaults to \code{FALSE}.}
}
\details{
The bagplot has been proposed by Rousseeuw et al. (1999) as a generalisation of the boxplot to bivariate data. It is constructed based on halfspace depth and as such is invariant under affine transformations. Similar graphical representations can be obtained by means of other depth functions, as illustrated in Hubert and Van der Veeken (2008) and in Hubert et al. (2015). See \code{\link{compBagplot}} for more details.

The deepest point is indicated with a "+" sign, the outlying observations with red stars.

The plot is made using \code{ggplot2}. The plot itself is returned by the function and is fully customisable using standard \code{ggplot2} commands.
}
\seealso{
  \code{\link{compBagplot}}, \code{\link{hdepth}}, \code{\link{projdepth}}, \code{\link{sprojdepth}}.
}
\references{
Rousseeuw P.J., Ruts I., Tukey, J.W. (1999). The bagplot: a bivariate boxplot. \emph{The American Statistician}, \bold{53}, 382--387.

Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.

Hubert M., Rousseeuw P.J., Segaert, P. (2015). Rejoinder to 'Multivariate functional outlier detection'. \emph{Statistical Methods & Applications}, \bold{24}, 269--277.

}

\author{P. Segaert}

\examples{
data(bloodfat)

# The bagplot can be plotted for the halfspace depth,
# the projection depth or the skewness-adjusted projection depth.
# Note that the projection depth is not appropiate for skewed data.
bagplot(compBagplot(bloodfat))
bagplot(compBagplot(bloodfat, type = "projdepth"))
bagplot(compBagplot(bloodfat, type = "sprojdepth"))

# The mean features of the bagplot can easily be adjusted.
result <- compBagplot(bloodfat)
bagplot(result, databag = FALSE, dataloop = FALSE)
bagplot(result, colorbag = rgb(0.2,0.2,0.2), colorloop = "green")


data(cardata90)
result <- compBagplot(cardata90)
bagplot(result)

# Compared to the original paper on the bagplot,
# an additional outlier is identified. However this
# point lies very close to the fence and this may be
# attributed to differences in numerical rounding.
# This may be illustrated by plotting the fence.
plot <- bagplot(result, plot.fence = TRUE)
plot

# The returned object is a ggplot2 object and may be
# edited using standard ggplot2 commands.
library("ggplot2")
plot + ylab("Engine displacement") + xlab("Weight in pounds")


}
\keyword{Graphical}
