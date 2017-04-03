\name{stars}
\alias{stars}
\title{stars data}
\description{

 Data corresponding to a Hertzsprung-Russel diagram of the star cluster CYG OB1 containing 47 stars in the direction of Cygnus.
 A typical Hertzsprung-Russel diagram shows the logarithm of the temperature in reverse order from high to low.

 The data has been taken from Humphreys (1978) by C. Doom who calibrated them according to Vansina en De Greve (1982).

}
\usage{
data(stars)
}
\format{
Following variables are included:
 \describe{
    \item{\code{LogTemp}}{Logarithm of the effective temperature at the star's surface.}
    \item{\code{LogLight}}{Logarithm of a star's light intensity. }
  }
}
\source{
Humphreys R.M. (1978). Studies of luminous stars in nearby galaxies. I. Supergiants and O stars in the milky way.  \emph{Astrophysics Journal Supplementary Series}, \bold{38}, 309--350.

Hand D.J., Daly F., Lunn A., McConway A. (1994). A Handbook of Small Data Sets. \emph{Londen: Chapman and Hall},  dataset 367.
}
\references{
Vansima F., De Greve J.P. (1982). Close binary systems before and after mass transfer. \emph{Astrophysics and Space Science}, \bold{87}, 377--401. 

Rousseeuw P.J., Leroy A. (1987). Robust Regression and Outlier Detection. \emph{New York: Wiley}. 
}
\examples{
data(stars)
plot(stars,xlim=rev(range(stars[,1])))
}
\keyword{datasets}
