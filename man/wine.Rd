\name{wine}
\alias{wine}
\docType{data}
\title{Proton Nuclear Magnetic Resonance spectra of 40 different wine samples}
\description{
The original data set consists of Proton Nuclear Magnetic Resonance (NMR) spectra of 40 different wine samples in the spectral region from 6.00 to 0.50. This data set corresponds to the region between wavelengths 5.62 and 5.37 only for which \eqn{t = 397} measurements are available for each curve. The data has been analyzed in Hubert et al. (2015), see below. 
}
\usage{data("wine")}
\format{
  A three dimensional \eqn{t  = 397} by \eqn{n = 40} by \eqn{p = 1} array,
                with \eqn{t} the number of observed time points,
                \eqn{n} the number of functional observations
                and \eqn{p} the number of measurements
                for every functional observation at every wavelength.
}
\details{
When using this data set, please cite both of the references below. 
}
\source{
Larsen F, van den Berg F, Engelsen S (2006) An exploratory chemometric study of NMR spectra of table
wines. \emph{Journal of Chemometrics}, \bold{20} - (5), 198-208}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection (with rejoinder). \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}
\examples{
data(wine)
matplot(wine[,,1], type="l", lty=1, col = "black")
}
\keyword{datasets}
