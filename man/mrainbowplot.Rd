\name{mrainbowplot}
\alias{mrainbowplot}
\title{Rainbow plot for bivariate data}
\description{
Makes a scatterplot of bivariate data and colors the observations according to their depth value.
}
\usage{
mrainbowplot(x, depths, col = NULL, plot.options = list())
}
\arguments{
  \item{x}{An \eqn{n} by 2 data matrix.
         }
  \item{depths}{A column vector of length \eqn{n}. \cr
                The depth values of the observations in \code{x}. The coloring
                is based on these depth values.}
  \item{col}{An \eqn{m > 2} by 3 matrix. \cr
             Colors in rgb format. The user may use this argument to set the colorscale of the depth range. The first row should contain the rgb values for the lowest depth value, the last row the rgb values of the color for the deepest depth value. Colors for other depth values are interpolated. When more then two rows are provided the color range will be equidistantly divided over the different colors.
             }
  \item{plot.options}{A list of available options:
                      \itemize{
                      \item{legend.title} Title of the legend. \cr
                        Defaults to \code{"Depth"}.
                      \item{point.size} Numeric defining the size of the points in the plot. \cr
                        Defaults to 4.
                      }
                     }
}
\details{
The plot is made using ggplot2. The plot itself is returned by the function and is fully customisable using standard ggplot2 commands. Similar plots for multivariate data with \eqn{p > 2} can be made using the ggpairs function in the library \code{GGally}. 
}

\author{P. Segaert}
\examples{
data(cardata90)
Result <- projdepth(x = cardata90)
plot.options <- list(legend.title = "PD")
plot <- mrainbowplot(cardata90,
                     depths = Result$depthZ, plot.options = plot.options)
library("ggplot2")
plot + ggtitle("Rainbowplot of the cardata using projection depth.")

#The default color range may be adjusted using the col argument.
RGBmatrix <- c(1, 0, 0, #Red
               1, 1, 1, #White
               0, 1, 0) #Green
RGBmatrix <- matrix(RGBmatrix, ncol = 3, byrow = TRUE)
plot <- mrainbowplot(cardata90, 
                     depths = Result$depthZ, col = RGBmatrix,
                     plot.options = plot.options)
plot + ggtitle("Rainbowplot of the cardata using projection depth.")

}

\keyword{Graphical}

