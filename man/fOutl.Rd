\name{fOutl}
\alias{fOutl}
\title{Functional outlyingness measures for functional data}
\usage{
  fOutl(x, z = NULL, type = "fAO", alpha = 0, time = NULL, 
        diagnostic = FALSE, distOptions = NULL)
}
\description{
Computes several measures of functional outlyingness for multivariate 
functional data. 
}

\arguments{
       \item{x}{A three dimensional \eqn{t} by \eqn{n} by \eqn{p} array,
                with \eqn{t} the number of observed time points,
                \eqn{n} the number of functional observations
                and \eqn{p} the number of measurements
                for every functional observation at every time point.}
       \item{z}{An optional three-dimensional \eqn{t} by \eqn{m} by \eqn{p} 
                array, containing the observations for which to compute the 
                functional outlyingness with respect to \code{x}. 
                If \code{z} is not specified, it is set equal to \code{x}. 
                The time points of \code{z} should correspond to those of 
                \code{x}. }
    \item{type}{The depth used in the computations.
                The outlyingness measure used in the computations.
                One of the following options: \code{"fSDO"},
                \code{"fAO"}, \code{"fDO"} or \code{"fbd"}. \cr
                Defaults to \code{"fAO"}.}
   \item{alpha}{Specifies the weights at every cross-section. 
                When \code{alpha = 0}, uniform weights are used. Otherwise
                \code{alpha} should be a weight vector of length \eqn{t}. \cr
                Defaults to 0. }
    \item{time}{If the measurements are not equidistant, 
                a sorted numeric vector containing a set of time points. \cr
                Defaults to \code{1:t}.}
\item{diagnostic}{If set to \code{TRUE}, the output contains some additional 
                  components: \cr 
                  crossDists: an \eqn{n} by \eqn{t} matrix containing the 
                  multivariate outlyingness of each observation at each time 
                  point \cr
                  locOutl: output containing flags for local outlyingness 
                  (see "Value" for more details)
                  \cr
                  Defaults to FALSE.
                  }
\item{distOptions}{A list of options to pass to the function
                   calculating the cross-sectional distances. \cr
                    See \code{outlyingness}, \code{adjOutl}, or
                    \code{bagdistance}.}
}
\details{
The functional outlyingness of a multivariate curve with respect to a 
given set of multivariate curves is defined as the weighted average of its 
multivariate outlyingness at each time point (Hubert et al., 2015). 
The functional outlyingness can be computed in all dimensions \eqn{p} using
bagdistance, projection depth and skewness-adjusted projection depth. 

When the data array \code{z} is specified, the functional outlyingness and diagnostic 
information for the data array \code{x} is also returned whenever the underlying
outlyingness routine allows it. For more information see the specific routines 
listed in the section "See Also"

In some situations, additional diagnostics are available to flag outlying time 
points. At each time point, observations from the data array \code{x} are marked
if they are flagged as outliers. The observations from the data array \code{x} 
are marked if their scaled outlyingness is larger than a prescribed cut-off value
from the chi-square distribution. For more details see the respective 
outlyingness routines

It is possible that at certain time points a part of the algorithm can not be 
executed due to e.g. exact fits. In that case the weight of that particular time
point is set to zero. A warning is issued at the end of the algorithm to signal 
these time points. Furthermore the output contains an extra argument giving the 
indices of the time points where problems occured.

}
\value{
 A  list with the following
  components:
            \item{fOutlyingnessX}{Vector of length \eqn{n} containing the 
                                  functional outlyingness of every curve
                                  from \code{x}.}
            \item{fOutlyingnessZ}{Vector of length \eqn{m} containing the 
                                  functional outlyingness of every curve
                                  from \code{z}.}
            \item{weights}{Vector of weights according to the input parameter 
                           \code{alpha}.}
 \item{crossDistsX}{An \eqn{n} by \eqn{t} matrix containing the multivariate
                    outlyingness of each observation of \code{x} at each point. 
                    Only provided if the input parameter \code{diagnostic} is 
                    set to \code{TRUE}.}
 \item{crossDistsZ}{An \eqn{m} by \eqn{t} matrix containing the multivariate
                    outlyingness of each observation of \code{z} at each point. 
                    Only provided if the input parameter \code{diagnostic} is 
                    set to \code{TRUE}.}
 \item{locOutlX}{An \eqn{n} by \eqn{t} matrix flagging local outlyingness for 
                \code{x}. Only provided if the input parameter \code{diagnostic}
                is set to TRUE.\cr
                The \eqn{(i,j)}th element takes value 1 if curve 
                \eqn{x_i} is outlying at time point \eqn{j}.
               }
 \item{locOutlZ}{An \eqn{m} by \eqn{t} matrix flagging local outlyingness for 
                \code{z}. Only provided if the input parameter \code{diagnostic}
                is set to TRUE.\cr
                The \eqn{(i,j)}th element takes value 1 if curve 
                \eqn{z_i} is outlying at time point \eqn{j}.
               }
 \item{IndFlagExactFit}{Vector containing the indices of the
                        time points for which an exact fit is
                        detected.}
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection.
\emph{Statistical Methods and Applications,} \bold{24}(2), 177--202.

Hubert M., Rousseeuw P.J., Segaert P. (2016). Multivariate and functional classification using depth and distance.
\emph{Advances in Data Analysis and Classification,} in press.

}
\seealso{\code{\link{bagdistance}}, \code{\link{outlyingness}}, \code{\link{adjOutl}}}
\examples{
# We will illustrate the function using a univariate functional sample.
data(octane)
Data <- octane

# When the option diagnostic is set to TRUE, a crude diagnostic
# to detect outliers can be extracted from the local outlyingness
# indicators. 
Result <- fOutl(x = Data, type = "fAO", diagnostic = TRUE)
matplot(Data[,,1], type = "l", col = "black", lty = 1)
for (i in 1:dim(Data)[2]) {
  if(sum(Result$locOutlZ) > 0) {
    obsData <- matrix(Data[,i,1], nrow = 1)
    obsData[!Result$locOutlZ[i,]] <- NA
    obsData <- rbind(obsData, obsData)
    matpoints(t(obsData), col = "red", pch = 15)
  }
}
# For more advance outlier detection techniques, see the 
# fom routine.
}
\author{P. Segaert}
\keyword{functional}


