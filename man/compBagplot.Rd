\name{compBagplot}
\alias{compBagplot}
\title{Location estimate based on the halfspace depth}
\description{
Computes all elements of the bagplot, a generalisation
of the univariate boxplot to bivariate data. The bagplot can be
computed based on halfspace depth, projection depth or skewness-adjusted projection depth. To make the actual plot, the function \code{bagplot} needs to be called on the result of \code{compBagplot}.
}
\usage{
compBagplot(x, type = "hdepth", sizesubset = 500,
            extra.directions = FALSE, options = NULL)
}
\arguments{
	\item{x}{An \eqn{n} by 2 data matrix.}

  \item{type}{Determines the depth function used to construct the bagplot: \code{"hdepth"} for halfspace depth, \code{"projdepth"} for projection depth and \code{"sprojdepth"} for skewness-adjusted projection depth. \cr
              Defaults to \code{"hdepth"}.
}
	\item{sizesubset}{When computing the bagplot based on halfspace depth,
                    the size of the subset used to perform the main
                    computations. See Details for more information. \cr
                    Defaults to \eqn{500}.
                    }
  \item{extra.directions}{Logical indicating whether additional directions should
                         be considered in the computation of the fence for the
                         bagplot based on projection depth or skewness-adjusted
                         projection depth. If set to \code{TRUE} an additional 
                         250 equispaced directions
                         are added to the directions defined by the points in
                         \code{x} themselves and the center. If \code{FALSE}
                         only directions determined by the points in \code{x} are
                         considered. \cr
                         Defaults to \code{FALSE}.
                        }
  \item{options}{A list of options to pass to the \code{projdepth} or \code{sprojdepth} function.
                 In addition the following option may be specified:
                  \itemize{
                      \item{ \code{max.iter} \cr
                            The maximum number of iterations in the bisection algorithm used
                            to compute the depth contour corresponding to the cutoff. See
                            \code{depthContour} for more information. \cr
                            Defaults to \eqn{100}.
                           }
                          }

                }

}
\details{
The bagplot has been proposed by Rousseeuw et al. (1999) as a generalisation of the boxplot to bivariate data.  It is constructed based on halfspace depth. In the original format the deepest point is indicated by a "+" and is contained in the bag which is defined as the depth region containing the 50\% observations with largest depth. The fence is obtained by inflating the bag (relative to the deepest point) by a factor of three. The loop is the convex hull of the observations of \code{x} inside the fence. Observations outside the fence are flagged as outliers and plotted with a red star. This function only computes all the components constituting the bagplot. The bagplot itself can be drawn using the \code{bagplot} function.

The bagplot may also be defined using other depth functions. When using projection depth or skewness-adjusted projection depth the bagplot is build as follows. The center corresponds to the observation with largest depth. The bag is constructed as the convex hull of the fifty percent points with largest depth. Outliers are identified as points with a depth smaller than a cutoff value, see \code{projdepth} and \code{sprojdepth} for the precise definition.
The loop is computed as the convex hull of the non-outlying points. The fence is approximated by the convex hull of those points that lie on rays from the center through the vertices of the bag and have a depth that equals the cutoff depth. For a better approximation the user can set the input parameter \code{extraDirections} to \code{TRUE} such that an additional 250 equally spaced directions on the circle are considered.

The computation of the bagplot based on halfspace depth can be time
consuming. Therefore it is possible to limit the bulk of the computations
to a random subset of the data. Computations of the halfspace median and
the bag are then based on this random subset. The number of points in this
subset can be controlled by the optional argument \code{sizesubset}.

It is first checked whether the data is found to lie on a line. If so, the routine will give a warning, giving back the dimension of the subspace (being 1) together with the normal vector to that line.


}
\value{
 A  list with components: \cr
	\item{center}{Center of the data. \cr
                       When \code{type = "hdepth"}, this corresponds with the
                       Tukey median.
                       In other cases this point corresponds to the point
                       with maximal depth.
	}
  \item{chull}{When \code{type = "hdepth"}, these are the vertices of
               the region with maximal halfspace depth.
               In other cases this is a null vector.
  }
  \item{bag}{The coordinates of the vertices of the bag.}
  \item{fence}{The coordinates of the vertices of the fence.}
  \item{datatype}{An \code{(nx3)} matrix. The first two columns correspond with \code{x}. The third column indicates the
                 position of each observation of \code{x} in the bagplot: 2 for observations in
                 the bag, 1 for the observations in the fence and 3 for
                 outliers. \cr
                 Note that points may not be in the same order as in \code{x}.

  }
	\item{flag}{A vector of length \code{n} wich is 0 for outliers and 1 for regular observations of \code{x}.}
  \item{depth}{The depth of the observations of \code{x}.}
  \item{dimension}{If the data are lying in a lower dimensional subspace, the dimension of this subspace.}
  \item{hyperplane}{If the data are lying in a lower dimensional subspace, a direction orthogonal to this subspace.}
  \item{type}{Same as the input parameter \code{type}.}
}
\seealso{
  \code{\link{bagplot}}, \code{\link{hdepth}}, \code{\link{projdepth}}, \code{\link{sprojdepth}}.
}
\references{
Rousseeuw P.J., Ruts I., Tukey J.W. (1999). The bagplot: A bivariate boxplot. \emph{The American Statistician}, \bold{53}, 382--387.

Hubert M., Van der Veeken S. (2008). Outlier detection for skewed data. \emph{Journal of Chemometrics}, \bold{22}, 235--246.

Hubert M., Rousseeuw P.J., Segaert, P. (2015). Rejoinder to 'Multivariate functional outlier detection'. \emph{Statistical Methods & Applications}, \bold{24}, 269--277.
}

\author{P. Segaert based on Fortran code by P.J. Rousseeuw, I. Ruts and A. Struyf.}

\examples{
data(bloodfat)
Result <- compBagplot(bloodfat)
bagplot(Result)

# The sizesubset argument may be used to control the
# computation time when computing the bagplot based on
# halfspace depth. However results may be unreliable when
# choosing a small subset for the main computations.
system.time(Result1 <- compBagplot(bloodfat))
system.time(Result2 <- compBagplot(bloodfat, sizesubset = 100))
bagplot(Result1)
bagplot(Result2)

# When using the projection depth or skewness-adjusted
# projection depth to compute the bagplot a list of options
# may be passed down to the (adj)outlyingness routines.
options <- list(type = "Rotation",
                ndir = 50,
                stand = "unimcd",
                h = floor(nrow(bloodfat)*3/4))
Result <- compBagplot(bloodfat,
                      type = "projdepth", options = options)
bagplot(Result)

# The fence is computed using the depthContour function.
# To get a smoother fence, one may opt to consider extra
# directions.
options <- list(ndir = 500,
                seed = 36)
Result <- compBagplot(bloodfat,
                      type = "sprojdepth", options = options)
bagplot(Result, plot.fence = TRUE)

options <- list(ndir = 500,
                seed = 36)
Result <- compBagplot(bloodfat,
                      type = "sprojdepth", options = options,
                      extra.directions = TRUE)
bagplot(Result, plot.fence = TRUE)
}
\keyword{multivariate}
