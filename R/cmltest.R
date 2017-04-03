cmltest <- function(x, z){

  ######
  # Check input.
  if (missing(x)) {
    stop("Input argument x is required.")
  }
  if (missing(z)) {
    stop("Input argument z is required.")
  }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {
    stop("The input argument x must be a numeric data matrix.")
  }
  n1 <- nrow(x)
  p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }
  if (p1 != 2) {
    stop("Only applicable to bivariate data.")
  }
  # Check the z data.
  z <- data.matrix(z)
  if (!is.numeric(z)) {
    stop("The input argument z must be a numeric data matrix.")
  }
  n2 <- nrow(z)
  p2 <- ncol(z)
  if (p1 != p2) {
    stop("Data dimension has to be the same for x and z.")
  }
  if (n2 > sum(complete.cases(z))) {
    stop("Missing values in z are not allowed.")
  }
  if (n2 != 1) {
    stop("z can only contain one observation.")
  }

  result <- rdepth(x = x, z = z)

  k <- floor(result$depthZ[1] * n1)
  if (k > round((n1 - 1) / 2)) {
    pvalue <- 1
  } else {
    jprime <- 0:round(k / (n1 - 2 * k))
    pvalue <- 2 * (n1 - 2 * k) * sum(dbinom(x = (n1 - k + jprime * (n1 - 2 * k)),
                                            size = n1,
                                            prob = 0.5)
    )
  }

  return(pvalue)
}
