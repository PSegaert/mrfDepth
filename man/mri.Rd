\name{mri}
\alias{mri}
\docType{data}
\title{Intensities of MRI images}
\description{
Felipe et al. (2005) obtained intensities of MRI images of 9 different parts of the human
body (plus a group consisting of all remaining body regions, which was of course very
heterogeneous). They then transformed their data to univariate curves. }
\usage{data("plane")}
\format{
A list of arrays corresponding to each bodypart. For each bodypart, a three-dimensional \eqn{t = 99} by \eqn{n} by \eqn{p = 1} array is available. The index \eqn{t} corresponds to the different points of measurement, the index \eqn{n} to the different observations. 
}
\details{
When using this data set please cite both Felipe et al. (2005) and Hubert et al. (2017).
}
\source{
Felipe J.C., Traina A.J.M., Traina C. (2005). Global warp metric distance: boosting
content-based image retrieval through histograms. Proceedings of the Seventh IEEE
International Symposium on Multimedia (ISM05), p.8.

Chen, Y., Keogh, E., Hu, B., Begum, N., Bagnall, A., Mueen, A., Batista, G.J. (2015). The
UCR Time Series Classification Archive. [http://www.cs.ucr.edu/~eamonn/time_series_data] 
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2017). Multivariate and functional classification using depth and distance. \emph{Advances in Data Analysis and Classification}, 11(3), 445-466.
}
\examples{
data(mri)
par(mfrow = c(2,1))
matplot(y = mri$bodypart1[,,1],
        type ="l",col = "black", lty = 1, xlab = "", ylab="x-coordinate", main = "plane 1")
matplot(y = mri$bodypart2[,,1],
        type ="l",col = "black", lty = 1, xlab = "", ylab="x-coordinate", main = "plane 2")
par(mfrow = c(1,1))
}
\keyword{datasets}
