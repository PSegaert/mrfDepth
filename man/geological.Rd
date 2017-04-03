\name{geological}
\alias{geological}
\title{Composition of elements in soil samples around the Baltic Sea}
\description{

This data originates from a geological survey on the composition in agricultural soils from 10 countries surrounding the Baltic Sea.
Top soil (0-25 cm) and bottom soil (50-75 cm) samples from 768 sites were analysed. This data frame contains the measurements corresponding to the total concentration of four elements in the top soil samples.
}
\usage{
data(geological)
}
\format{
A dataframe conraining following variables:
 \describe{
      \item{\code{Fe2O3}}{ Iron Oxide }
      \item{\code{MgO}}{ Magnesium oxide }
      \item{\code{MnO}}{ Manganese oxide }
      \item{\code{TiO2}}{ Titanium dioxide }
 }
}
\source{
Reimann C., Siewers U., Tarvainen T., Bityukova L., Eriksson J., Gilucis A.,  Gregorauskiene V., Lukashev V., Matinian N.N., Pasieczna A. (2000). Baltic soil survey: total concentrations of major and selected trace elements in arable soils from 10 countries around the Baltic Sea. \emph{Science of the Total Environment},
\bold{257}, 155--170.
}
\examples{
data(geological)
plot(geological)
}
\keyword{datasets}
