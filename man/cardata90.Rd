\name{cardata90}
\alias{cardata90}
\title{Car data from Consumer Reports in 1990}
\description{
  Subset from data on cars taken from pages 235-255, 281-285 and 287-288 of the April 1990 Consumer Reports Magazine.
}
\usage{
data(cardata90)
}
\format{
 A dataframe conraining following variables:
    \describe{
      \item{\code{Weight}}{ Weight of the car. }
      \item{\code{Disp}}{ Engine displacement in cubic inches. }
    }
}
\source{
\emph{Consumer Reports}, April 1990, 235--288.

Chambers J.M., Hastie T.J. (1993). Statistical Models in S. \emph{Londen: Chapman and Hall}, 46--47.
}
\references{
Rousseeuw, P.J., Ruts I, and Tukey J.W. (1999).
The bagplot: a bivariate boxplot.
\emph{The American Statistician} \bold{53}(4), 382--387.
}
\examples{
data(cardata90)
plot(cardata90)
}
\keyword{datasets}
