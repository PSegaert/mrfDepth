\name{mfd}
\alias{mfd}
\title{Multivariate functional depth for functional data}
\usage{
  mfd(x, z = NULL, type = "hdepth", alpha = 0, time = NULL, diagnostic = FALSE,
     depthOptions = NULL)
}
\description{
Computes the multivariate functional depth for multivariate functional data.
}

\arguments{
       \item{x}{A three dimensional \eqn{t} by \eqn{n} by \eqn{p} array,
                with \eqn{t} the number of observed time points,
                \eqn{n} the number of functional observations
                and \eqn{p} the number of measurements
                for every functional observation at every time point.}
       \item{z}{An optional three-dimensional \eqn{t} by \eqn{m} by \eqn{p} 
                array, containing the observations for which to compute the 
                multivariate functional depth with respect to \code{x}. 
                If \code{z} is not specified, it is set equal to \code{x}. 
                The time points of \code{z} should correspond to those of 
                \code{x}. }
    \item{type}{The depth used in the computations.
                One of the following options: \code{"hdepth"},
                \code{"projdepth"}, \code{"sprojdepth"}, \code{"sdepth"}. \cr
                Defaults to \code{"hdepth"}.}
   \item{alpha}{Specifies the weights at every cross-section. 
                When \code{alpha = 0}, uniform weights are used. Weights 
                following equation (2) in Claeskens et al. (2014) are obtained 
                by setting \code{alpha} to a number smaller than the maximal 
                depth at any time point. The weights are then proportional to 
                the volume of the \eqn{\alpha}-depth regions at each 
                cross-section. Otherwise \code{alpha} should be a weight vector 
                of length \eqn{t}. \cr
                Defaults to 0. }
    \item{time}{If the measurements are not equidistant, 
                a sorted numeric vector containing a set of time points. \cr
                Defaults to \code{1:t}.}
\item{diagnostic}{If set to \code{TRUE}, the output contains some additional 
                  components: \cr 
                  crossDepths: an \eqn{n} by \eqn{t} matrix containing the 
                  multivariate depth of each observation at each time point \cr
                  locOutl: output containing flags for local outlyingness 
                  (see "Value" for more details)
                  \cr
                  Defaults to FALSE.}
\item{depthOptions}{A list of options to pass to the function
                    calculating the cross-sectional depths. \cr
                    See \code{hdepth}, \code{projdepth}, \code{sprojdepth},
                    or \code{sdepth}.}
}
\details{
The multivariate functional depth of a multivariate curve with respect to a 
given set of multivariate curves is defined as the weighted average of its 
multivariate depth at each time point (Claeskens et al., 2014). 
The MFD can be computed in all dimensions \eqn{p} using halfspace depth, 
projection depth and skewness-adjusted projection depth. 
For \eqn{p<=2} also simplicial depth is available. 

When the data array \code{z} is specified, the MFD depth and diagnostic 
information for the data array \code{x} is also returned whenever the underlying
depth routine allows it. For more information see the specific depth routines 
listed in the section "See Also"

For the weight vector, three options are available: uniform weights, 
user-defined weights or weights proportional to the volume of the 
\eqn{\alpha}-depth contour at each time point. 
The \eqn{\alpha}-depth contours are computed using the \code{depthContour} function. 

In some situations, additional diagnostics are available to flag outlying time 
points, as described in Hubert et al. (2012). At each time point, observations 
from the data array \code{x} are marked if they are flagged as outliers. 
When using any of the projection depth measures, this flag is automatically 
returned by the corresponding functions. When using halfspace depth, 
the diagnostic is only available for bivariate curves. The observations 
from the data array \code{x} are marked if they are flagged as outliers by the 
bagplot, or similarly if their \code{bagdistance} is larger than 3 at that time 
point. This can be seen as a measure of local outlyingness. The option is not 
available for simplicial depth.  

It is possible that at certain time points a part of the algorithm can not be 
executed due to e.g. exact fits. In that case the weight of that particular time
point is set to zero. A warning is issued at the end of the algorithm to signal 
these time points. Furthermore the output contains an extra argument giving the 
indices of the time points where problems occured.


}
\value{
 A  list with the following  components:
 \item{MFDdepthX}{Vector of length \eqn{n} containing the MFD depth
                  of every curve from \code{x}.}
 \item{MFDdepthZ}{Vector of length \eqn{m} containing the MFD depth
                  of every curve from \code{z}.}
   \item{weights}{Vector of weights according to the input parameter 
                  \code{alpha}.}
 \item{crossDepthsX}{An \eqn{n} by \eqn{t} matrix containing the multivariate
                    depth of each observation of \code{x} at each point. Only 
                    provided if the input parameter \code{diagnostic} is set to 
                    \code{TRUE}.}
 \item{crossDepthsZ}{An \eqn{m} by \eqn{t} matrix containing the multivariate
                    depth of each observation of \code{z} at each point. Only 
                    provided if the input parameter \code{diagnostic} is set to 
                    \code{TRUE}.}
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
 \item{IndFlagBag}{Vector containing the indices of the
                        time points for which the bagplot could
                        not be computed.}
 \item{IndFlagIso}{Vector containing the indices of the time
                    points for which the cross-sectional \eqn{\alpha}-depth 
                    contours could not be computed.}
 \item{IndFlagAlpha}{Vector containing the indices of the
                      time points for which the volume of the
                      cross-sectional \eqn{\alpha}-depth contours could not
                      be computed.}
}
\references{
Claeskens G., Hubert M., Slaets L., Vakili K. (2014). Multivariate functional halfspace depth.
\emph{Journal of the American Statistical Association,} \bold{109}, 505, 411--423.
 
Hubert M., Claeskens G., De Ketelaere B., Vakili K. (2012).
A new depth-based approach for detecting outlying curves. In ,
In \emph{Proceedings of COMPSTAT 2012, edited by A. Colubi, K. Fokianos, G. Gonzalez-Rodriguez, E.J. Kontoghiorghes,} 329--340. 
 
}

\seealso{\code{\link{depthContour}}, \code{\link{hdepth}}, \code{\link{projdepth}}, \code{\link{sprojdepth}}, \code{\link{sdepth}}}

\examples{
data(octane)
Result <- mfd(x = octane, alpha = 0.125, diagnostic = TRUE)

}
\author{P. Segaert and M. Hubert}
\keyword{functional}


