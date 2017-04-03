\name{hdepth}
\alias{hdepth}
\title{Halfspace depth of points relative to a dataset}
\description{
Computes the halfspace depth of \eqn{p}-dimensional points \code{z} relative to a \eqn{p}-dimensional dataset \code{x}. Computation is exact for \eqn{p \le 3} and approximate when \eqn{p > 3}. For the approximate algorithm  the halfspace depth is computed as the minimal univariate halfspace depth over many directions. To obtain the univariate halfspace depth in the direction \eqn{v}, the dataset \code{x} is projected on \eqn{v}, and the univariate location depth of the points of \eqn{v'z_i} to \eqn{xv} is computed.
}
\usage{
hdepth(x, z = NULL, options = list())
}
\arguments{
  \item{x}{An \eqn{n} by \eqn{p} data matrix with observations in the rows and variables in the columns.}
  \item{z}{An optional \eqn{m} by \eqn{p} matrix containing
           rowwise the points \eqn{z_i} for which to compute
           the halfspace depth. If \code{z} is not specified,
           it is set equal to \code{x}.}
  \item{options}{ A list of available options:
                  \itemize{
                      \item{ \code{type} \cr
                      Determines the desired type of invariance for the approximate algorithm and should be
                      one of \code{"Affine"}, \code{"Rotation"} or \code{"Shift"}.
                      When the option \code{"Affine"} is used,  the directions
                      \eqn{v} are orthogonal to hyperplanes spanned by \eqn{p}
                      observations from \code{x}. When the option \code{"Rotation"}
                      is used, the directions pass by two randomly selected
                      observations from \code{x}. With the option \code{"Shift"},
                      directions are randomly generated.\cr
                      Defaults to \code{"Affine"}.}
                      \item{ \code{ndir} \cr
                      Determines the number of directions \eqn{v} by setting
                      \code{ndir} to a specific number or to \code{"all"}.
                      In the latter case, an exhaustive search over all possible
                      directions (according to \code{type}) is performed.
                      If \code{ndir} is larger than the number of possible directions,
                      the algorithm will automatically use this setting. \cr
                      Defaults to
                      \eqn{250p} when \code{type="Affine"}, to 5000 when
                      \code{type="Rotation"} and to 12500 when \code{type="Shift"}.}
                      \item{ \code{approx} \cr
                      The user may force approximate calculation in two or three
                      dimensions by setting this option to \code{TRUE}. \cr
                      Defaults to \code{FALSE}.}
                      \item{ \code{seed} } \cr
                      A strictly positive integer specifying the seed to be used
                      by the \code{C++} code.\cr
                      Defaults to \eqn{10}.
                  }
                }
}
\details{
Halfspace depth has been introduced by Tukey (1975). The halfspace depth of a point \eqn{z_i} is defined as the minimal number of observations from \code{x} that are contained in any closed halfspace with boundary through \eqn{z_i}.

In dimensions \eqn{p=2} and \eqn{p=3} the computations are by default carried out exactly using the algorithms described in Rousseeuw and Ruts (1996) and Rousseeuw and Struyf (1998). This yields an affine invariant measure of depth.
Approximate algorithms are also implemented which depending on the value chosen for \code{type} are affine, rotation or shift invariant. They can be used in any dimension. The shift invariant algorithm coincides with the random Tukey depth (Cuesta-Albertos and Nieto-Reyes, 2008).

The resulting halfspace depth values are invariant to affine transformations when the exact algorithm is used and invariant to affine transformations, rotations and shifts depending on the choice for \code{type} provided that the \code{seed} is kept fixed at different runs of the algorithm. Note that the halfspace depth values values are guaranteed to increase when more directions are considered provided the seed is kept fixed, as this ensures that the random directions are generated in a fixed order.

If the halfspace depth needs to be computed for \eqn{m} points \eqn{z_i}, it is recommended to apply the function once with the matrix \eqn{z} as input, instead of applying it \eqn{m} times with input vectors \eqn{z_i}, as numerous computations can be saved. The approximate algorithms automatically then also compute the depth values of the observations in \code{x}.
When only the halfspace depth of the observations in \code{x} is required, the call to the function should be 'hdepth(x)' or equivalently 'hdepth(x,x)'. In that case the depth values will be stored in the 'depthZ' output field. For bivariate data these will be the exact values by default.
To visualize the depth of bivariate data one can apply the \code{mrainbowplot} function. It plots the data colored according to their depth.

It is first checked whether the data lie in a subspace of dimension smaller than \eqn{p}. If so, a warning is given, as well as the dimension of the subspace and a direction which is orthogonal to it.
}
\value{
 A  list with components: \cr
      \item{depthX}{Vector of length \eqn{n} giving the halfspace depth of the observations in \code{x}.\cr
      By default exact if \eqn{p \le 3} and approximate if \eqn{p > 3} or the option \code{approx} is set to \code{TRUE}.}
  \item{depthZ}{Vector of length \eqn{m} giving the halfspace depth of the points in \code{z} relative to \code{x}. \cr
  By default exact if \eqn{p \le 3} and approximate if \eqn{p > 3} or the option \code{approx} is set to \code{TRUE}.
        		            }
      \item{singularSubsets}{When the input parameter type is equal to \code{"Affine"},
                             the number of \eqn{p}-subsets that span a subspace of
                             dimension smaller than \eqn{p-1}. In that case the orthogonal
                             direction can not be uniquely determined. This is an indication
                             that the data are not in general position.
                             When the input parameter type is equal to \code{"Rotation"} it
                             is possible that two randomly selected points of the data coincide
                             due to ties in the data. In this case this value signals how many
                             times this happens.
                             }
      \item{dimension}{When the data \code{x} are lying in a lower dimensional subspace,
                       the dimension of this subspace.
                       }
      \item{hyperplane}{When the data \code{x} are lying in a lower dimensional subspace,
                        a direction orthogonal to this subspace.
                        }
}
\seealso{
  \code{\link{hdepthmedian}}, \code{\link{mrainbowplot}}, \code{\link{bagdistance}}, \code{\link{bagplot}}
}
\references{
Tukey (1975). Mathematics and the picturing of data. \emph{Proceedings of the International Congress of Mathematicians,} \bold{2}, 523--531, Vancouver.

Rousseeuw, P.J., Ruts, I. (1996). AS 307: Bivariate location depth. \emph{Journal of the Royal Statistical Society: Series C},  \bold{45}, 516--526.

Rousseeuw, P.J., Struyf, A. (1998). Computing location depth and regression depth in higher dimensions. \emph{Statistics and  Computing}, \bold{8}, 193--203.

Cuesta-Albertos, J. and Nieto-Reyes, A. (2008). The random Tukey depth. \emph{Computational Statistics & Data Analysis}, \bold{52}:4979--4988.
}
\author{P. Segaert based on Fortran code by P.J. Rousseeuw, I. Ruts and A. Struyf, and \code{C++} code by P. Segaert and K. Vakili.}
\examples{
# Compute the halfspace depth of a simple
# two-dimensional dataset. 
data(cardata90)
Result <- hdepth(x = cardata90)
mrainbowplot(cardata90, depths = Result$depthZ)

# In two dimensions we may also opt to use the
# approximate algorithm. The number of directions
# may be specified through the option list.
options <- list(type = "Rotation",
                ndir = 750,
                approx = TRUE)
Result <- hdepth(x = cardata90, options = options)
# The resulting halfspace depth is monotone decreasing 
# in the number of directions.
options <- list(type = "Rotation",
                ndir = 10,
                approx = TRUE)
Result1 <- hdepth(x = cardata90, options = options)
options <- list(type = "Rotation",
                ndir = 500,
                approx = TRUE)
Result2 <- hdepth(x = cardata90, options = options)
which(Result1$depthZ - Result2$depthZ < 0)
# This is however not the case when the seed is changed
options <- list(type = "Rotation",
                ndir = 10,
                approx = TRUE)
Result1 <- hdepth(x = cardata90, options = options)
options <- list(type = "Rotation",
                ndir = 50,
                approx = TRUE,
                seed = 897)
Result2 <- hdepth(x = cardata90, options = options)
which(Result1$depthZ - Result2$depthZ < 0)
plot(Result1$depthZ - Result2$depthZ,
     xlab = "Index", ylab = "Difference in halfspace depth")

# We can also consider directions through two data
# points. If the sample is small enough one may opt
# to search over all choose(n,2) directions.
# Note that the computational load increases dramatically
# as n becomes larger.
options <- list(type = "Rotation",
                ndir = "all",
                approx = TRUE)
Result1 <- hdepth(x = cardata90, options = options)

# Alternatively one may consider randomly generated directions.
options <- list(type = "Shift",
                ndir = 250,
                approx = TRUE)
Result1 <- hdepth(x = cardata90, options = options)
}
\keyword{multivariate}
