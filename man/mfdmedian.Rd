\name{mfdmedian}
\alias{mfdmedian}
\title{Multivariate functional median for functional data.}
\description{
Computes the multivariate functional median, an estimate for the central tedency
for multivariate functional data.
}
\usage{
  mfdmedian(x, type="hdepth", crossDepthsX = NULL,
            depthOptions = NULL, centerOption = "maxdepth")
}
\arguments{
       \item{x}{A three dimensional \eqn{t} by \eqn{n} by \eqn{p} array,
                with \eqn{t} the number of observed time points,
                \eqn{n} the number of functional observations
                and \eqn{p} the number of measurements
                for every functional observation at every time point.}
    \item{type}{The depth used in the computations.
                One of the following options: \code{"hdepth"},
                \code{"projdepth"}, \code{"sprojdepth"}, \code{"dprojdepth"},
                \code{"sdepth"}. \cr  Defaults to \code{"hdepth"}.}
\item{crossDepthsX}{Depths at each time point. Can be used to save time}
\item{depthOptions}{A list of options to pass to the function
                    calculating the cross-sectional depths. \cr
                    See \code{hdepth}, \code{projdepth}, \code{sprojdepth},
                    \code{dprojdepth} or \code{sdepth}.}
\item{centerOption}{By default the 'maxdepth' center estimate is used for the
                     projmedian and the sprojmedian. This can be changed by setting
                     centerOption to 'gravity' or 'Huber' depending on wheter
                     projmedian or sprojmedian is used.}
}
\details{
The multivariate functional median a multivariate functional data set is defined
as the multivariate curve connecting the cross-sectional multivariate depth 
medians (Claeskens et al., 2014). 
The MFD median can be computed in all dimensions \eqn{p} using halfspace depth, 
projection depth and skewness-adjusted projection depth. 

It is possible that at certain time points a part of the algorithm can not be 
executed due to e.g. exact fits. In that case the estimate for the center will 
be set to NaN. A warning is issued at the end of the algorithm to signal 
these time points. Furthermore the output contains an extra argument giving the 
indices of the time points where problems occured.
}
\value{
 A  list with the following
  component:
       \item{MFDmedian}{An \code{(txp)} matrix containing the estimated central curve.}
       \item{IndFlagExactFit}{Vector containing the indices of the
                              time points for which an exact fit is
                              detected.}
}
\references{
Claeskens G., Hubert M., Slaets L., Vakili K. (2014). Multivariate functional halfspace depth.
\emph{Journal of the American Statistical Association,} \bold{109}, 505, 411--423.
}
\examples{
# We will illustrate the function using a bivariate functional sample.
data(characterTrajectories)
Data <- characterTrajectories[,1:50,]
Result <- mfdmedian(Data)

par(mfrow = c(1,2))
matplot(Data[,,1], type = "l", col = "black", lty = 1, ylab = "x-coordinate")
matlines(Result$MFDmedian[,1], type = "l", col = "red", lty = 1)
matpoints(Result$MFDmedian[,1], col = "red", pch = 15)
matplot(Data[,,2], type = "l", col = "black", lty = 1, ylab = "y-coordinate")
matlines(Result$MFDmedian[,2], type = "l", col = "red", lty = 1)
matpoints(Result$MFDmedian[,2], col = "red", pch = 15)
par(mfrow = c(1,1))

# Other depth function such as adjusted outlyingness may also 
# used to determine the cross-sectional depth median. 
# In this case the depth median is calculated by the 
# sprojmedian routine. 
# When different depth median are available, a specific depth
# median can be selected usinf the centerOption. Optional 
# arguments used by the sprojmedian routine may be specified
# using the depthOptions. For example one might choose the
# "Rotation" equivariance for 300 directions. 
Result <- mfdmedian(Data, type = "sprojdepth",
                    depthOptions = list(type = "Rotation",
                                        ndir = 300),
                    centerOption = "gravity")
par(mfrow = c(1,2))
matplot(Data[,,1], type = "l", col = "black", lty = 1, ylab = "x-coordinate")
matlines(Result$MFDmedian[,1], type = "l", col = "red", lty = 1)
matpoints(Result$MFDmedian[,1], col = "red", pch = 15)
matplot(Data[,,2], type = "l", col = "black", lty = 1, ylab = "y-coordinate")
matlines(Result$MFDmedian[,2], type = "l", col = "red", lty = 1)
matpoints(Result$MFDmedian[,2], col = "red", pch = 15)
par(mfrow = c(1,1))

# If the user already placed a call to the mdf routine
# with the diagnostic options set to TRUE, the 
# mfdmedian can easily be obtained by passing the cross-sectional 
# depths. This drastically saves computing time.  
tResult <- mfd(x = Data, type = "sprojdepth", diagnostic = TRUE)
Result <- mfdmedian(Data, type = "sprojdepth",
                    crossDepthsX = tResult$crossdepthX,
                    centerOption = "gravity")
}
\author{P. Segaert}
\keyword{functional}
