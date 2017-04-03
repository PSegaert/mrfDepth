hdepth <- function(x, z = NULL, options = list()){

  ######
  # Check input.
  if (missing(x)) {
    stop("Input argument x is required.")
  }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {
    stop("The input argument x must be a numeric data matrix.")
  }
  n1 <- nrow(x)
  p <- p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }
  # Check the z data.
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
    stop("The data dimension has to be the same for x and z.")
  }
  if (n2 > sum(complete.cases(z))) {
    stop("Missing values in z are not allowed.")
  }

  # check options en load defaults if necessary.
  if (is.null(options)) {
    options <- list()
  }
  if (!is.list(options)) {
    stop("options must be a list")
  }
  if ("ndir"  %in% names(options)) {
    ndir <- options[["ndir"]]
  } else {
    ndir <- NULL
  }
  if ("type" %in% names(options)) {
    type <- options[["type"]]
  } else {
    type <- "Affine"
  }
  if ("approx" %in% names(options)) {
    approx <- options[["approx"]]
  } else {
    approx <- FALSE
  }
  if ("seed"  %in% names(options)) {
    seed <- options[["seed"]]
  } else {
    seed <- NULL
  }


  # Check number of directions and type.
  type.id <- match(type, c("Affine", "Rotation", "Shift"))[1]
  if (is.na(type.id)) {
    stop("The input parameter type must be one of: Affine, Rotation or Shift.")
  }
  if (is.null(ndir)) {
    if (type.id == 1)  ndir <- 250 * ncol(x)
    if (type.id == 2)  ndir <- 250 * 20
    if (type.id == 3)  ndir <- 250 * 50
  }
  calc.exact <- 0
  # If the specified number of directions is greater than
  # the possible different directions, switch to exact calculation.
  if (is.numeric(ndir)) {
    if (ndir < 1) {
      stop("The number of directions must be a positive integer.")
    }
    if (type.id == 1) {
      ndir0 <- choose(n1, p1)
      if (ndir0 <= ndir) {
       ndir <- ndir0
       calc.exact <- 1
      }
    }
    if (type.id == 2) {
      ndir0 <- choose(n1, 2)
      if (ndir0 <= ndir) {
        ndir <- ndir0
        calc.exact <- 1
      }
    }
  }
  if (!is.numeric(ndir)) {
    if (ndir == "all") {
      if (type.id == 1) {
        ndir <- choose(n1, p1)
        calc.exact <- 1
        if (ndir > 1e7) {
          stop("ndir is larger than 1e7. Try a smaller value of ndir.")
        }
      }
      if (type.id == 2) {
        ndir <- choose(n1, 2)
        calc.exact <- 1
        if (ndir > 1e7) {
          stop("ndir is larger than 1e7. Try a smaller value of ndir.")
        }
      }
      if (type.id == 3) {
        stop("Cannot compute all directions for type Shift.")
      }
    } else {
      stop("The input parameter ndir is not recognized.")
    }
  }
  if (n1 < (p1 + 1) & type.id == 1) {
    stop("When type is affine, n should be larger than p.")
  }
  if (is.null(seed)) {
    seed <- 10
  }
  if (!is.numeric(seed)) {
    stop("The seed must be a strictly positive integer.")
  }
  if (seed <= 0) {
    stop("The seed must be a strictly positive integer.")
  }

  #####
  # Check data for possible exact fit situations.
  tol <- 1e-7
  scaled.x <- scale(x)
  temp <- attributes(scaled.x)
  column.sd <- temp[["scaled:scale"]]
  if (sum(column.sd <= 1e-14) > 0) {
    warning("One of the variables has zero standard deviation.
            Check the data matrix x.")
    returned.result <- list(depthX = NULL,
                            depthZ = NULL,
                            singularSubsets = NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14))
    class(returned.result) <- c("mrfDepth", "hdepth")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n1 - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x
            could not be computed.")
    returned.result <- list(depthX = NULL,
                            depthZ = NULL,
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "hdepth")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more information.")
    returned.result <- list(depthX = NULL,
                            depthZ = NULL,
                            singularSubsets = NULL,
                            dimension = sum(w1$d > tol),
                            hyperplane = w1$v[, which(w1$d == min(w1$d))[1]])
    class(returned.result) <- c("mrfDepth", "hdepth")
    return(returned.result)
  }

  #####
  #Do the calculations

  else if (p == 2 & approx == FALSE) {
    result <- .Fortran("HSDEP2",
                     as.double(z[, 1, drop = TRUE]), #1 First coordinate
                                                     #  of points to compute
                                                     #  the depth of.
                     as.double(z[, 2, drop = TRUE]), #2 Second coordinate of
                                                     #  points to compute the
                                                     #  depth of.
                     as.integer(n2),                 #3 Number of points in
                                                     #  the set z.
                     as.double(x[, 1, drop = TRUE]), #4 First coordinate of the
                                                     #  data matrix x.
                     as.double(x[, 2, drop = TRUE]), #5 Second coordinate of the
                                                     #  data matrix x.
                     as.integer(n1),                 #6 Number of observations
                                                     #  in x.
                     as.double(rep(0, n2)),          #7 Halfspace depth of the
                                                     #  points in z.
                     as.double(rep(0, n2)),          #8 Simplicial depth of the
                                                     #  points in z.
                    PACKAGE = "mrfDepth")
    returned.result <- list(depthX = NULL,
                            depthZ = result[[7]],
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
  }
  else if (p == 3 & approx == FALSE) {
    result <- .Fortran("HSDEP3",
                     as.double(z[, 1, drop = TRUE]), #1 First coordinate of
                                                     #  points to compute the
                                                     #  depth of.
                     as.double(z[, 2, drop = TRUE]), #2 Second coordinate of
                                                     #  points to compute the
                                                     #  depth of.
                     as.double(z[, 3, drop = TRUE]), #3 Third coordinate of
                                                     #  points to compute the
                                                     #  depth of.
                     as.integer(n2),                 #4 Number of points in
                                                     #  the set z.
                     as.double(x[, 1, drop = TRUE]), #5 First coordinate of
                                                     #  the data matrix x.
                     as.double(x[, 2, drop = TRUE]), #6 Second coordinate of
                                                     #  the data matrix x.
                     as.double(x[, 3, drop = TRUE]), #7 Third coordinate of
                                                     #  the data matrix x.
                     as.integer(n1),                 #8 Number of observations
                                                     #  in x.
                     as.double(rep(0, n2)),          #9 Halfspace depth of the
                                                     #  points in z.
                     as.double(rep(0, n2)),          #10 Flag signaling exact
                                                     #   fit.
                     PACKAGE = "mrfDepth")
    returned.result <- list(depthX = NULL,
                            depthZ = result[[9]],
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
  }
  else{
    result <- .C("HSDND",
               as.integer(n1 + n2),         #1 Total number of points for x
                                            #  and z combined.
               as.integer(p),               #2 Dimension of the data.
               as.integer(ndir),            #3 Number of directions.
               as.double(rbind(x, z)),      #4 Data matrix (both x and z).
               as.integer(rep(1, n1 + n2)), #5 Halfspace depth of all points
               as.integer(0),               #6 Number of singular directions.
               as.integer(type.id),         #7 Integer indicating which type
                                            #  of directions to consider.
               as.integer(n1),              #8 Number of points in the data
                                            #  matrix x.
               as.integer(calc.exact),      #9 Flag indicating whether exact
                                            #  computation should be performed.
               as.integer(seed),            #10 The random seed.
               PACKAGE = "mrfDepth")
    returned.result <- list(depthX = (result[[5]][1:n1] / n1),
                            depthZ = (result[[5]][(n1 + 1):(n1 + n2)]) / n1,
                            singularSubsets = result[[6]],
                            dimension = NULL,
                            hyperplane = NULL)
  }
  class(returned.result) <- c("mrfDepth", "hdepth")
  return(returned.result)
}
