\name{fom}
\alias{fom}
\title{Draws the Functional Outlier Map (FOM)}
\description{
Creates the Functional Outlier Map, a graphical tool to detect outliers in multivariate functional data. Depending on the position of a multivariate curve in the FOM, different types of outliers may be distinguished. 
}
\usage{
fom(fOutlResult, cutoff = FALSE)
}
\arguments{
  \item{fOutlResult}{The return of a call to \code{fOutl}}
  \item{cutoff}{Boolean indicating whether the cutoff should be drawn.}
}
\details{
The fom is only applicable when \code{fOutl} is called with the following options: \code{diagnostic = TRUE}.

The functional outlier map was proposed by Hubert et al. (2015) and subsequently improved by Rousseeuw et al. (2017). It consists of graphical tool to detect outliers in multivariate functional data based on functional outlyingness measures such as the \eqn{fAO} or \eqn{fDO} (see \code{fOutl}). Similar grpahical displays may be obtained by considering other functional outlyingness measures such as the \eqn{fBD} or \eqn{fSDO}. 

The coordinates of the points in the FOM correspond to two outlyingness indicators. On the horizontal axis, the functional outlyingness as obtained by the routine \code{\link{fOutl}} is plotted. On the vertical axis the scaled standard deviation of the cross-sectional outlyingness measures across time are plotted. The FOM thus consists of the following points representing the curve \eqn{i}: \eqn{(fOutl_i; std_j(fOutl_i(t_j)) / (1 + fOutl_i))}. The scaling of the standard deviation is added to ensure that curves with a different location, measured from the center, but with the same relative variability have a similar \eqn{y}-coordinate. 

For some underlying multivariate outlyingness measures, the return of \code{fOutl} contains an additional argument signalling whether curve \eqn{i} is outlying as measured in the multivariate space at time point \eqn{j}. This information is incorperated in the FOM by plotting the points with a different degree of outlyingness with a different symbol. Curves that are outlying in the multivariate space determined by at least \eqn{75\%} of the total number of time points in the domain are plotted using filled diamonds. Similarly curves that are outlying in at least \eqn{50\%} or \eqn{25\%} of he time points are plotted in filled squares and filled triangles respectively. Curves that are outlyingn in fewer than \eqn{25\%} of the time points plotted using filled circles. 

The user may opt to draw a cut off line for the detection of outliers on the FOM.
This cutoff was introduced in Rousseeuw et al. (2017) and is based on the euclidean
distance between the points on the FOM and the origin, after scaling with the median.
More specifically, let \eqn{(fOutl_i;vOutl_i) = (fOutl_i; std_j(fOutl_i(t_j)) / (1 + fOutl_i))} and 
let \eqn{cOutl_i^2 = (fOutl_i / median(fOutl))^2 + (vOutl_i / median(vOutl))^2}.
Finally, with \eqn{lcOutl_i = \log(0.1 + coutl_i)},an observation lies outside of the cutoff when 
\eqn{(lcOutl_i - median(lcOutl)) / mad(lcOutl) > \Phi^{-1}(0.995)}.

This FOM may be read in a way similar to the outlier map of robust regression (Rousseeuw and van Zomeren 1990) and the outlier map of
robust principal components (Hubert et al. 2005). Points in the lower left part of the FOM represent regular curves which hold a central position in the dataset. Points in the lower right part are curves with a high \eqn{fOutl} but a low variability of their cross-sectional \eqn{fOutl} values. This happens for shift outliers, i.e. curves which have the same shape as the majority but are shifted on the whole domain. Points in the upper left part have a low \eqn{fOutl} but a high variability in cross-sectional \eqn{fOutl}. Typical examples are isolated outliers, i.e. cuves which only display outlyingness over a small part of their domain.
The points in the upper right part of the FOM have both a high \eqn{fOutl} and a high cross-sectional \eqn{fOutl}. These correspond to curves which are strongly outlying on a substantial part of their domain. 
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection (with rejoinder). \emph{Statistical Methods & Applications}, \bold{24}, 177--202.

Rousseeuw P.J., Raymaekers J., Hubert M. (2017). A Measure of Directional Outlyingness with Applications to Image Data and Video. \emph{arXiv}, arXiv:1608.05012v2.
}

\author{P. Segaert}

\examples{
data(octane)

# To construct the FOM, one first need to calculate 
# the functional outlyingness.
# Note that the option diagnostic in fOutl must be 
# set to TRUE. If not calling fom will result in an 
# error
Result <- fOutl(octane, alpha = 0, type = "fAO", diagnostic = TRUE);
fom(Result)

# The user may opt to draw a cut off line seperating the outliers.
# which will be plotted in red
fom(Result, cutoff = TRUE)

# Six observations are flagged as outliers. These correspond to 
# the sample with added ethanol. 

}
\keyword{Graphical}

