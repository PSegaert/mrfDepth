\name{distSpace}
\alias{distSpace}
\title{distSpace}
\description{
Calculation of distance space representation.
}
\usage{
distSpace(trainingData, testData = NULL, type = "bagdistance", options = NULL)
}
\arguments{
  \item{trainingData}{A list of \code{nxp} matrices containing the multivariate data or a list of \eqn{t} by \eqn{n} by \eqn{p} arrays containing the functional data.}
  \item{testData}{An \code{mxp} matrix containing all multivariate training data or a \eqn{t} by \eqn{m} by \eqn{p} array for functional data.}
  \item{type}{The distance used in the computations.
                For multivariate data one of the following options: \code{"bagdistance"},
                \code{"outlyingness"}, or \code{"adjOutl"}. \cr
                For functional data one of the following options: \code{"fBD"}, \code{"fSDO"} or \code{"fAO"}
                Defaults to \code{"bagdistance"}.}
  \item{options}{ A list of options to pass to the function
                  calculating the underlying distance. \cr
                  See \code{"bagdistance"}, \code{"outlyingness"}, \code{"adjOutl"} or \code{fOutl} for more information.}
}
\details{
The distance is a tool in supervised classification and was introduced in Hubert et al. (2016) as a generalisation of the depth-depth representation of a multivariate sample. Based on a distance transform, an observation (be it multivariate or functional) is mapped to its representation in distance space. The distance transformation consists of mapping the observation to a vector containing at coordinate \eqn{i} the distance to the training group \eqn{i}. After transformation, any multivariate classifier may be used to classify new observations in distance space. Typically the \eqn{k}-nearest neighbour algorithm is used. 

Different options are available to calculate the distance to each of the training groups. For multivariate data, the user may choose between the bagdistance or any of the projection type distance including the Stahel-Donoho outlyingness, the adjusted outlyingness or the directional outlyingness. For functional data, the user may opt to employ the functional bagdistance (fBD), the functional Stahel-Donoho (fSDO)the functional adjusted-outlyingness (fAO) or the functional directional outlyingness (fDO). Options to available in each of the underlying distance routines may be passed down using the \code{options} argument.
}
\value{
A matrix \code{q x (p+1)} composed of two blocks. The first block contains in each row an observations in the training set with in each column the distance to each of the training sets. The last column contains a label indicating the original group membership of the observation. The second block contains the observations in the test set, if any, with in each column the distance to the different training sets. The last column now contains an inddicator signalling the observations was part of the test set. 
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2017). Multivariate and functional classification using depth and distance. \emph{Advances in Data Analysis and Classification}, 11(3), 445-466.
}
\author{P. Segaert}

\examples{

# We will use two multivariate toy data sets
data(cardata90)
data(bloodfat)

# Build the training data
trainingData <- list(set1 = cardata90,
                     set2 = bloodfat)
# Transform the data into distspace
Result <- distSpace(trainingData = trainingData)
# Plot the results
plotColors <- c(rep("orange", nrow(cardata90)),
                rep("blue", nrow(bloodfat)))
plot(Result[, 1:2],
     col = plotColors,
     xlab = "distance to cardata90", ylab = "distance to bloodfat",
     main = "distspace representation of cardata90 and the bloodfat data.")


# By default the bagdistance is used to transform the data. 
# This can be changed by using the type argument. Additional option to be
# passed to the underlying function calculatin the distance may be passed in 
# the option argument.
options <- list(type = "Affine", ndir = 1000, seed = 3)
Result <- distSpace(trainingData = trainingData,
                    type = "adjOutl",
                    options = options)
# Plot the results
plotColors <- c(rep("orange", nrow(cardata90)),
                rep("blue", nrow(bloodfat)))
plot(Result[, 1:2],
     col = plotColors,
     xlab = "distance to cardata90", ylab = "distance to bloodfat",
     main = "distspace representation of cardata90 and the bloodfat data.")

data(octane)
data(glass)
trainingData <- list(set1 = glass[1:100,, , drop = FALSE],
                     set2 = octane[1:100,, , drop = FALSE])
# Transform the data into distspace
Result <- distSpace(trainingData = trainingData, type = "fAO")

}

\keyword{functional}
\keyword{multivariate}
