\name{characterI}
\alias{characterI}
\title{Writing trajectories of the letter i.}
\description{
Subset of the 'Character Trajectories Data Set' from the UCI Machine Learning Repository. The data set consists of trajectories of the tip of a pen whilst writing the letter 'i'. All samples are from the same writer. Original data has been processed.
}
\usage{
data(characterI)
}
\format{
Three dimensional array. The first dimension represents time. The second dimension corresponds to the observation number. The third dimensions contain the X and Y coordinates.
}
\source{
Bache K. and Lichman, M. (2013). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science.
}
\references{
Williams B.H., Toussaint M. and Storkey A.J. (2006) Extracting motion primitives from natural handwriting data. In ICANN, volume 2, pages 634--643, .

Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection (with rejoinder). \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}
\examples{
data(characterI)
par(mfrow = c(1,2))
matplot(y = characterI[,,1],
        type ="l",col = "black", lty = 1, xlab = "Time", ylab="X position of the pen")
matplot(y = characterI[,,2],
        type ="l",col = "black", lty = 1, xlab = "Time", ylab="Y position of the pen")
par(mfrow = c(1,1))
}
\keyword{datasets}
