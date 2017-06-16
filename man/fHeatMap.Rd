\name{fheatmap}
\alias{fheatmap}
\title{Creates MFOD heatmap}
\description{
Creates the heatmap in paper either for depth or for distances.
}
\usage{
fheatmap(rowValues,cellValues,type,scalename="")
}
\arguments{
  \item{rowValues}{a}
  \item{cellValues}{a}
  \item{type}{ One of \code{"depth"} or \code{"distance"}. Determines whether a depth or a distance heatmap is made.}
  \item{scalename}{a}
}
\details{fg
}
\references{
Hubert M., Rousseeuw P.J., Segaert P. (2015). Multivariate functional outlier detection (with rejoinder). \emph{Statistical Methods & Applications}, \bold{24}, 177--202.
}

\author{P. Segaert}

\examples{
library("ggplot2")
data(octane)
Result <- mfd(octane, diagnostic = TRUE, type ="sprojdepth")
Plot <- fheatmap(rowValues = Result$MFDdepthZ,
                 cellValues = Result$crossdepthZ,
                 type = "depth",
                 scalename ="SPD")
                 
Result <- fOutl(octane, diagnostic = TRUE, type ="fAO")
Plot <- fheatmap(rowValues = Result$fOutlyingnessZ,
                 cellValues = Result$crossDistsZ,
                 type = "distance",
                 scalename ="AO")                 
                 
#Customize the look of the plot
#Plot <- Plot + xlab("time point")
#Plot
}
\keyword{Graphical}
