\name{bloodfat}
\alias{bloodfat}
\title{Blood data for patients with narrowing arteries}
\description{

Data were collected on the concentration of plasma cholesterol and plasma triglycerides (mg/dl) for 371 male patients evaluated for chest pain. For 51 of those patients, no evidence of heart disease was found. This subset corresponds to the remaining 320 patients for which there was evidence of narrowing arteries.

}
\usage{
data(bloodfat)
}
\format{
A dataframe conraining following variables:
 \describe{
      \item{\code{Cholesterol}}{Concentration of plasma cholesterol [mg/dl].}
      \item{\code{Triglycerides}}{Concentration of plasma triglycerides [mg/dl].}
      
 }
}
\source{
Hand D.J., Daly F., Lunn A., McConway A. (1994). A Handbook of Small Data Sets. \emph{Londen: Chapman and Hall}, dataset 277.
}
\references{
Scott D.W., Gotto A.M., Cole J.S., Gorry G.A. (1978). Plasma lipids as collateral risk factors in coronary artery disease: a study of 371 males with chest pain.  \emph{Journal of Chronic Diseases}, \bold{31}, 337--345.
}
\examples{
data(bloodfat)
plot(bloodfat)
}
\keyword{datasets}
