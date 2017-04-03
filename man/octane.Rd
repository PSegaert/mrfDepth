\name{octane}
\alias{octane}
\title{Near infrared spectra of gasoline samples.}
\description{
The last six observations contain added ethanol.
}
\usage{
data(octane)
}
\format{
Matrix observations correspond to rows.
}
\source{
Esbensen K. (2001). Multivariate data analysis in practice. \emph{5th edn. Camo Software}, Trondheim, Norway
}
\references{
Rousseeuw P.J., Debruyne M., Engelen S., Hubert M. (2006). Robustness and outlier detection in chemometrics. \emph{
Critical Reviews in Analytical Chemistry}, \bold{36}, 221--242.
}
\examples{
data(octane)
matplot(octane[,,1], type="l", lty=1, col = "black")
}
\keyword{datasets}
