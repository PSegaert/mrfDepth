\name{plotContours}
\alias{plotContours}
\title{Draws depth contours of bivariate data}
\description{
Draws the depth contours of bivariate data computed with \code{depthContour}.

}
\usage{
plotContours(x, depthContour, data = TRUE)
}
\arguments{
\item{x}{An \eqn{n} by 2 data matrix.}
\item{depthContour}{ The result of a call to \code{depthContour}. }
\item{data}{Logical value indicating whether the data \code{x} should be plotted.}
}
\details{
The plot is made using \code{ggplot2}. The plot itself is returned by the function and is fully customisable using standard \code{ggplot2} commands.
}
\references{
Ruts I., Rousseeuw P.J. (1996). Computing depth contours of bivariate point clouds. \emph{Computational Statistics & Data Analysis}, \bold{23}, 153-168.

}
\seealso{
 \code{\link{depthContour}} 
}
\author{P. Segaert}

\examples{
data(cardata90)
Result <- depthContour(x=cardata90, alpha=c(0.02,0.125,0.3) , type="hdepth")
plotContours(x = cardata90, depthContour = Result)
Result <- depthContour(x=cardata90, alpha=c(0.1, 0.2, 0.3, 0.4) , type="projdepth")
plotContours(x = cardata90, depthContour = Result)
Result <- depthContour(x=cardata90, alpha=c(0.1, 0.2, 0.3, 0.4) , type="sprojdepth")
plotContours(x = cardata90, depthContour = Result)

}

\keyword{Graphical}

