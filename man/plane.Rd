\name{plane}
\alias{plane}
\docType{data}
\title{Fighter plane dataset}
\description{
The fighter plane dataset of Thakoor and Gao (2005) describes 7 shapes: of the Mirage,
Eurofighter, F-14 with wings closed, F-14 with wings opened, Harrier, F-22 and F-15.
Each class contains 30 shape samples obtained from digital pictures, which Thakoor and
Gao (2005) then reduced to the univariate functions. 
}
\usage{data("plane")}
\format{
A list of arrays corresponding to each plane. For each plane, a three-dimensional \eqn{t} by \eqn{n} by \eqn{p = 1} array is available. The index \eqn{t} corresponds to the different points of measurement, the index \eqn{n} to the different observation. 
}
\details{
When using this data set please cite both Thakoor et al. (2005) and Hubert et al. (2017).
}
\source{
Thakoor, N. and Gao, J. (2005). Shape classifier based on generalized probabilistic descent
method with hidden Markov descriptor. Tenth IEEE International Conference on
Computer Vision (ICCV 2005), Vol. 1: 495-502.

Chen, Y., Keogh, E., Hu, B., Begum, N., Bagnall, A., Mueen, A., Batista, G.J. (2015). The
UCR Time Series Classification Archive. [http://www.cs.ucr.edu/~eamonn/time_series_data] 
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2017). Multivariate and functional classification using depth and distance. \emph{Advances in Data Analysis and Classification}, 11(3), 445-466.
}
\examples{
data(plane)
par(mfrow = c(2,1))
matplot(y = plane$plane1[,,1],
        type ="l",col = "black", lty = 1, xlab = "", ylab="x-coordinate", main = "plane 1")
matplot(y = plane$plane2[,,1],
        type ="l",col = "black", lty = 1, xlab = "", ylab="x-coordinate", main = "plane 2")
par(mfrow = c(1,1))
}
\keyword{datasets}
