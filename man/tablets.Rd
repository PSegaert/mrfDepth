\name{tablets}
\alias{tablets}
\docType{data}
\title{Near Infrared Spectroscopy responses for a batch of pills}
\description{
The original data set consists of Near Infrared Spectroscopy (NIR) data of a batch of pills with different levels of the active compound.  This data set considers 70 oberservation weighing 90 mg and 20 observations weighing 250 mg and also differ in the amount of active compound. 

The data corresponds to the form discussed in Hubert et al. (2015), see below. 
}
\usage{data("wine")}
\format{
  A three dimensional \eqn{t  = 404} by \eqn{n = 90} by \eqn{p = 3} array,
                with \eqn{t} the number of observed time points,
                \eqn{n} the number of functional observations
                and \eqn{p} the number of measurements
                for every functional observation at every wavelength.
}
\details{
For each wavelength and each curve, a three-dimensional vector of observations is available. The second coordinate corresponds to the original data, the first coordinate is corresponds to the mean of the corresponding curve. The last component corresponds to the derivative of the original curve data.

When using this data, please cite both of the references listed below. 
}
\source{
Dyrby M., Engelsen S.B., Norgaard L. , Bruhn M.,  and Nielsen L. (2002).
Chemometric Quantitation of the Active Substance in a Pharmaceutical Tablet Using Near Infrared (NIR) Transmittance and NIR FT Raman Spectra
\emph{Applied Spectroscopy},  \bold{56} - (5).
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection (with rejoinder). \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}
\examples{
data(tablets)
par(mfrow=c(3,1))
matplot(tablets[,,1], type="l", lty=1, col = "black")
matplot(tablets[,,2], type="l", lty=1, col = "black")
matplot(tablets[,,3], type="l", lty=1, col = "black")
par(mfrow=c(1,1))
}
\keyword{datasets}
