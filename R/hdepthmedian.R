hdepthmedian <- function(x, maxdir = NULL) {

    ######
  # Check input.
  if (missing(x)) {
    stop("The input argument x is required.")
  }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {
    stop("The input argument x must be a numeric data matrix.")
  }
  n <- nrow(x)
  p <- ncol(x)
  if (n > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }

  #maxdir
  if (is.null(maxdir)) {
    maxdir <- 250 * p
  }
  if (!is.numeric(maxdir)) {
    stop("maxdir should be a numeric.")
  }
  if (maxdir < 1) {
    stop("maxdir should be at least 1.")
  }
  if (n < p + 1) {
    stop("There should a least be p+1 observations.")
  }

  # initialisation defaults
  nstp <- 5 * n ^ (0.3)
  ntry <- 10 * (p + 1)
  nalt <- 4 * (p + 1)

  #####
  # Check data for possible exact fit situations.
  # Check data for possible exact fit situations.
  tol <- 1e-7
  scaled.x <- scale(x)
  temp <- attributes(scaled.x)
  column.sd <- temp[["scaled:scale"]]
  if (sum(column.sd <= 1e-14) > 0) {
    warning("One of the variables has zero standard deviation.
            Check the data matrix x.")
    returned.result <- list(median = NULL,
                            depth = NULL,
                            dithered = NULL,
                            ndir = NULL,
                            AlgStopFlag = NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14))
    class(returned.result) <- c("mrfDepth", "hdepthmedian")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x
            could not be computed.")
    returned.result <- list(median = NULL,
                            depth = NULL,
                            dithered = NULL,
                            ndir = NULL,
                            AlgStopFlag = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "hdepthmedian")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more information.")
    returned.result <- list(list(median = NULL,
                                 depth = NULL,
                                 dithered = NULL,
                                 ndir = NULL,
                                 AlgStopFlag = NULL,
                                 dimension = sum(w1$d > tol),
                                 hyperplane = w1$v[, which.min(w1$d)]))
    class(returned.result) <- c("mrfDepth", "hdepthmedian")
    return(returned.result)
  }

  #####
  # Do the calculations

  if (p == 1) {
    returned.result <- list(median = median(x),
                            depth = min(length(which(x >= median(x))),
                                        length(which(x <= median(x)))) / n,
                            dithered = NULL,
                            ndir = NULL,
                            AlgStopFlag = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    return(returned.result)
  }
  else if (p == 2) {
    Result <- .Fortran("HALFMED2D",
              as.double(x[, 1, drop = TRUE]),  #1 First variable data matrix.
              as.double(x[, 2, drop = TRUE]),  #2 Second variable data matrix.
              as.integer(n),                   #3 Number of observations.
              as.double(rep(0, 2)),            #4 Halfspace median
              as.double(0),                    #5 Depth of halfspace median.
              as.integer(0),                   #6 Logical indicating dithering.
              PACKAGE = "mrfDepth")
    returned.result <- list(median = Result[[4]],
                            depth = Result[[5]],
                            dithered = as.logical(Result[[6]]),
                            ndir = NULL,
                            AlgStopFlag = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    return(returned.result)
  }
  else{
    Result <- .Fortran("HSDEPTH_DEEPEST",
                as.double(x),         #1 Data matrix.
                as.integer(n),        #2 Number of observations.
                as.integer(p),        #3 Number of dimensions.
                as.integer(maxdir),   #4 Maximum number of directions.
                as.integer(nstp),     #5 Maximum number of steps.
                as.integer(ntry),     #6 Maximum number of trials.
                as.integer(nalt),     #7 Maximum number of steps to increase
                                      #  depth.
                as.double(rep(0, p)), #8 Coordinates of deepest point.
                as.double(0),         #9 Halfspace depth of deepest point.
                as.integer(0),        #10 Logical indicating stopping criterium.
                as.integer(0),        #11 Number of directions used.
                as.integer(0),        #12 Number of steps used.
                PACKAGE = "mrfDepth")
    returned.result <- list(median = Result[[8]],
                            depth = Result[[9]],
                            dithered = NULL,
                            ndir  =  Result[[11]],
                            AlgStopFlag = Result[[10]],
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "hdepthmedian")
    return(returned.result)
  }

}
