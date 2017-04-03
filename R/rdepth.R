rdepth <- function(x, z = NULL, ndir = NULL){
  ######
  # Check input.
  if (missing(x)) {
    stop("Input argument x is required.")
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
  if (n1 < p1) {
    stop("At least p observations are required.")
  }
  #Check the z data.
  if (is.null(z)) {
    z <- x
  }
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
  # Older code requires intercept to be in the last column.
  z <- z[, c(2:(p2),1), drop = FALSE]
  
  #check ndir
  if (is.null(ndir)) {
    ndir <- 250 * p1
  }
  if (ndir < 1) {
    stop("The number of directions must be a positive integer.")
  }


  #####
  #Check data for possible exact fit situations.
  tol <- 1e-7
  scaled.x <- scale(x)
  temp <- attributes(scaled.x)
  column.sd <- temp[["scaled:scale"]]
  if (sum(column.sd <= 1e-14) > 0) {
    warning("One of the variables has zero
            standard deviation. Check the data matrix x.")
    returned.result <- list(depthZ = NULL,
                            singularSubsets = NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14)
                            )
    class(returned.result) <- c("mrfDepth", "rdepth")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n1 - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the
            data matrix x could not be computed.")
    returned.result <- list(depthZ = NULL,
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "rdepth")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit is found. Check the output for more details.")
    returned.result <- list(depthZ = NULL,
                            singularSubsets = NULL,
                            dimension = sum(w1$d > tol),
                            hyperplane = w1$v[, which(w1$d == min(w1$d))[1]]
                            )
    class(returned.result) <- c("mrfDepth", "rdepth")
    return(returned.result)
  }


  if (p1 == 2) {
    depth <- rep(NA, n2)
    for (i in 1:n2) {
      depth[i] <- rdepth2(z[i, ], x = x[, 1], y = x[, 2])
    }
    returned.result <- list(depthZ = depth / n1,
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "rdepth")
    return(returned.result)
  }
  else if (p1 == 3) {
    Fresult <- .Fortran("RDEPTH3",
                        as.double(z[, 1, drop = TRUE]), #1 slope of first
                                                        #  variable of Z
                        as.double(z[, 2, drop = TRUE]), #2 slope of second
                                                        #  variable of Z
                        as.double(z[, 3, drop = TRUE]), #3 intercepts of Z
                        as.integer(n2),                 #4 Number of planes in z
                        as.double(x[, 1, drop = TRUE]), #5 First variable of X
                        as.double(x[, 2, drop = TRUE]), #6 Second variable of X
                        as.double(x[, 3, drop = TRUE]), #7 Respons variable
                        as.integer(n1),                 #8 Number of planes in X
                        as.double(rep(0, n2)),          #9 Regression depths
                        as.integer(rep(0, n2)),         #10 Flag for exact fit.
                        PACKAGE = "mrfDepth")
    returned.result <- list(depthZ = Fresult[[9]],
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "rdepth")
    return(returned.result)
  }
  else if (p1 == 4) {
    Fresult <- .Fortran("RDEPTH4",
                        as.double(z),                   #1 Z variable
                        as.integer(n2),                 #2 Number of planes in z
                        as.double(x[, 1, drop = TRUE]), #3 First variable of X
                        as.double(x[, 2, drop = TRUE]), #4 Second variable of X
                        as.double(x[, 3, drop = TRUE]), #5 Third variable of X
                        as.double(x[, 4, drop = TRUE]), #6 Response variable
                        as.integer(n1),                 #7 Number of planes in X
                        as.double(rep(0, n2)),          #8 Regression depths
                        as.integer(rep(0, n2)),         #9 Flag for exact fit.
                        PACKAGE = "mrfDepth")
    returned.result <- list(depthZ = Fresult[[9]],
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "rdepth")
    return(returned.result)
  }
  else{
    storage.mode(x) <- "double"
    storage.mode(z) <- "double"
    Fresult <- .Fortran("RDEPTHND",
                        as.double(z),           #1 Z variable
                        as.integer(n2),         #2 Number of planes in z
                        as.double(x),           #3 X variables
                        as.integer(n1),         #4 Number of planes in X
                        as.integer(p1),         #5 Number of variables in X
                        as.integer(ndir),       #6 Number of directions to
                                                #  calculate
                        as.double(rep(0, n2)),  #7 Regression depths
                        as.integer(rep(0, n2)), #8 Number of singular directions
                        as.integer(rep(0, n2)), #9 Flag for all directions
                                                #  singular
                        PACKAGE = "mrfDepth")
    returned.result <- list(depthZ = Fresult[[7]],
                            singularSubsets = Fresult[[8]],
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "rdepth")
    return(returned.result)
  }
}

rdepth2 <- function(b, x, y, ordered = FALSE){
  # first value of b is the slope
  # sevond value of b is the intercept

  x <- as.vector(x)
  y <- as.vector(y)
  xy <- cbind(x, y)
  n <- nrow(xy)
  
  # sort the xy values
  if (!ordered) {
    xy <- xy[order(xy[, 1]), , drop = F]
    x <- xy[, 1]
    y <- xy[, 2]
  }

  resid <- y - b[1] * x - b[2]
  resid[abs(resid) < 10 ^ -7] <- 0
  
  positive.resid <- resid >= 0
  neg.resid <- resid <= 0
  
  dupx <- duplicated(x)
  
  if (sum(dupx)) {
    r1 <- (1:n)[!dupx]
    if (length(r1) == 1)
      r2 <- n - 1
    else r2 <- c(diff(r1) - 1, n - max(r1))
    r1 <- cbind(r1, r1 + r2)
    sumident <- function(x, y) {
      z <- apply(y[x[1]:x[2], , drop = F], 2, sum)
      return(z)
    }
    resid <- apply(r1, 1, sumident, cbind(positive.resid, neg.resid))
    positive.resid <- resid[1,  ]
    neg.resid <- resid[2,  ]
    n <- length(positive.resid)
  }
  lplus <- cumsum(positive.resid)
  rplus <- lplus[n] - lplus
  lmin <- cumsum(neg.resid)
  rmin <- lmin[n] - lmin
  depth <- pmin(lplus + rmin, rplus + lmin)
  return(min(depth))
}
