\encoding{UTF-8}
\name{cmltest}
\alias{cmltest}
\title{Test for linearity of the conditional median in simple regression}
\description{
A test based on regression depth for linearity of the
 conditional median \code{z} given the simple regression dataset \code{x}. The calculation of the regression depth of \code{z} with respect to \code{x} is done by the function \code{rdepth}.
 The test is only valid when \code{x} contains no duplicates.

}
\usage{
cmltest(x, z)
}
\arguments{
	\item{x}{An \eqn{n} by 2 regression data matrix. \cr
                    The first column is the explanatory variable, the second column corresponds to the response variable. }
  \item{z}{A parameter vector (intercept and slope). }
}
\details{
The following hypothesis test is performed:\cr
    \eqn{H0}: The data come from a model with: \eqn{med(x_2|x_1) = z_1 + z_2 * x_1} \cr
 The test statistic being used is the regression depth of \code{z} with respect to \code{x}.
}
\value{
	\item{pval}{The \eqn{p}-value of the hypothesis test.}
}
\seealso{
  \code{\link{rdepth}}, \code{\link{rdepthmedian}}
}
\references{
Van Aelst, Stefan, Rousseeuw, P.J., Hubert, M., Struyf, A. (2002). The Deepest Regression Method. \emph{Journal of Multivariate Analysis}, \bold{81}, 138--166.

Rousseeuw, P.J., Struyf, A. (2002). A Depth Test for Symmetry, in:  \emph{Goodness-of-Fit Tests and Model Validity}, \enc{Birkhäuser}{Birkhauser} Boston, 401--412.
}

\author{P. Segaert}

\examples{

data("stars")

# Compute the least squares fit. Due to outliers
# this fit will be bad and thus H0 should be rejected. 
temp <- lsfit(x = stars[,1], y = stars[,2])$coefficients
intercept <- temp[1]
slope <- temp[2]
z <- matrix(c(intercept, slope), nrow = 1)
pvalue1 <- cmltest(x = stars, z = z)
pvalue1

# Lets now test the deepest regression line. 
result <- rdepthmedian(x = stars)
pvalue2 <- cmltest(x = stars, z = matrix(result$deepest, nrow = 1))
pvalue2

plot(stars)
abline(a = intercept, b = slope)
abline(result$deepest, col = "red")
text(x = 3.8, y = 5.3, labels = paste("p-value", round(pvalue1, digits = 3)) )
text(x = 4.45, y = 4.8, labels = paste("p-value", round(pvalue2, digits = 3)), col = "red")


}
\keyword{regression}

