symtest <- function(x, z, options = list()){

  ######
  # Check input.
  if (missing(x)) {
    stop("Input argument x is required.")
  }
  if (missing(z)) {
    stop("Input argument z is required.")
  }

  #Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {
    stop("The input argument x must be a numeric data matrix.")
  }
  n1 <- nrow(x)
  p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }
  if (n1 < (p1 + 1)) {
    stop("At least (p+1) observations are needed.")
  }
  if (p1 != 2) {
    stop("Only applicable to bivariate data.")
  }
  

  #Check z argument
  if (!is.numeric(z)) {
    stop("The input argument z must be a numeric data matrix.")
  }
  z <- data.matrix(z)
  n2 <- nrow(z)
  p2 <- ncol(z)
  if (p2 == 1) {
    z <- t(z)
    n2 <- nrow(z)
    p2 <- ncol(z)
  }
  if (p1 != p2) {
    stop("Data dimension has to be the same for x and z.")
  }
  if (n2 > sum(complete.cases(z))) {
    stop("Missing values in z are not allowed.")
  }
  if (n2 != 1) {
    stop("Data matrix z must have one row.")
  }

  #check options
  if (is.null(options)) {
    options <- list()
  }
  if (!is.list(options)) {
    stop("options must be a list")
  }

  #Perform the test
  Temp <- hdepth(x = x, z = z, options = options)
  if (is.nan(Temp$depthZ[1])) {
    stop("The depth of z could not be calculated.
         Check the return of hdepth(x,z = z).")
  }

  k <- floor(Temp$depthZ[1] * n1)
  if (k > round((n1 - 1) / 2)) {
    pvalue <- 1
  }
  else {
      jprime <- 0:round(k / (n1 - 2 * k))
      pvalue <- 2 * (n1 - 2 * k) * sum(dbinom(x = (n1 - k + jprime * (n1 - 2 * k)),
                                              size = n1,
                                              prob = 0.5)
                                      )
  }
  return(pvalue)
}
