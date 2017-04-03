\name{glass}
\alias{glass}
\docType{data}
\title{EXPMA spectra of glass samples}
\description{
  The glass data set studied by Lemberge et al. (2000) consisting of 180 different 16-17th century archaeological glass samples.Electron Probe X-ray Microanalysis (EPXMA) intensities across \eqn{750} different wavelengths are recorded using a Jeol JSM 6300 scanning electron microscope equipped with an energy-dispersive Si(Li) X-ray detection system (SEM-EDX).
}
\usage{data("glass")}
\format{
  A three dimensional \eqn{t = 750} by \eqn{n = 180} by \eqn{p = 1} array.
}
\source{
  Lemberge, P., De Raedt, I., Janssens, K.H., Wei, F., and Van Espen, P.J. (2000).
  Quantitative Z-analysis of 16th-17th century archaeological glass vessels using
  PLS regression of EPXMA and \eqn{\mu}-XRF data. \emph{Journal of Chemometrics}, \bold{14},
  751--763.
}
\references{
  Hubert M., Rousseeuw P. J., and Vanden Branden K. (2005).
  ROBPCA: A New Approach to Robust Principal Component Analysis,
  \emph{Technometrics}, \bold{47}, 64--79.
}
\examples{
  data(glass)
  matplot(glass[,,1], type="l", lty=1, col = "black")
}