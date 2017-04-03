rdepthmedian <- function(x, maxit=NULL) {

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
  n1 <- nrow(x)
  p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }
  if (n1 < p1) {
    stop("At least p observations are required.")
  }
  if (p1 < 2) {
    stop("The input argument x should have at least two columns.")
  }
  #check maxit
  if (is.null(maxit)) {
    maxit <- 100
  }
  if (maxit < 1) {
    stop("The maximum number of iterations must be a positive integer.")
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
    returned.result <- list(deepest = NULL,
                            depth = NULL,
                            ndir = NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14)
                            )
    class(returned.result) <- c("mrfDepth", "rdepthmedian")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n1 - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the
            data matrix x could not be computed.")
    returned.result <- list(deepest = NULL,
                            depth = NULL,
                            ndir = NULL,
                            dimension = NULL,
                            hyperplane = NULL
                            )
    class(returned.result) <- c("mrfDepth", "rdepthmedian")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check output for more details.")
    returned.result <- list(deepest = NULL,
                            depth = NULL,
                            ndir = NULL,
                            dimension = sum(w1$d > tol),
                            hyperplane = w1$v[, which(w1$d == min(w1$d))[1]]
                            )
    class(returned.result) <- c("mrfDepth", "rdepthmedian")
    return(returned.result)
  }

  if (p1 == 2) {
    # Exact algorithm
    # All lines trough two points are considered.
    # Their depth is computed.
    # The deepest fit is the average of all fits with maximal depth.
    max.depth <- 0
    deepest.fit <- rep(0.0, 2)
    n.deepest <- 0

    x <- x[order(x[, 1]), , drop = F]

    for (i in 1:(n1 - 1)) {
      # Construct all lines through current data point and all data points
      # with a higher index. This way we avoid calculate the depth of lines
      # twice.
      numerator <- (x[(i + 1):n1, 2] - x[i, 2])
      denominator <- (x[(i + 1):n1, 1] - x[i, 1])
      slopes <- numerator / denominator
      intercepts <- (-1) * slopes * x[i, 1] + x[i, 2]
      cons.lines <- cbind(slopes, intercepts)
      # Lines corresponding to ties or verticals should not be considered.
      # If all considered lines should not be considered skip to next
      # iteration.
      rem.ind <- which(numerator == 0)
      rem.ind <- c(rem.ind, which(denominator == 0))
      if (length(rem.ind) > 0) {
        if (length(rem.ind) == nrow(cons.lines)) next()
        cons.lines <- cons.lines[-rem.ind, ]
      }
      # Calculate depth of all considered lines.
      depths <- rep(0.0, nrow(cons.lines))
      for (j in 1:nrow(cons.lines)) {
        temp.res <- rdepth2(b = cons.lines[j, ],
                            x = x[(i + 1):n1, 1], y = x[(i + 1):n1, 2],
                            ordered = TRUE)
        depths[j] <- temp.res
      }
      # Check if new maximal depth is found. If so reset the maximal depth and
      # deepest fit. Next check which lines have depth equal to maximum depth.
      # Add them to the sum in deepest fit and update counter on how many lines
      # were summed over.
      max.obtained <- max(depths, na.rm = TRUE)
      if (max.obtained > max.depth) {
        max.depth <- max.obtained
        deepest.fit <- rep(0.0, 2)
        n.deepest <- 0
      }
      if (max.obtained == max.depth) {
        ind.max <- which(depths == max.depth)
        deepest.fit <- deepest.fit +
                           colSums(cons.lines[ind.max, , drop = FALSE])
        n.deepest <- n.deepest + length(ind.max)
      }
    }
    # The deepest fit should be average, so divide by total number of lines in
    # the sommation.
    deepest.fit <- (1 / n.deepest) * deepest.fit
    deepest.fit <- rev(deepest.fit)
    names(deepest.fit) <- c("intercept", "slope")

    returned.result <- list(deepest = deepest.fit,
                            depth = max.depth / n1,
                            ndir = NULL,
                            dimension = NULL,
                            hyperplane = NULL
                            )
    class(returned.result) <- c("mrfDepth", "rdepthmedian")
    return(returned.result)

  } else {
    result <- .Fortran("SWEEPMEDRES",
                       as.double(x),          #1 The x data
                       as.integer(n1),        #2 Number of planes in x
                       as.integer(p1),        #3 Number of dimensions
                       as.double(rep(0, p1)), #4 Vector containing deepest
                       #  regression estimate
                       as.integer(maxit),     #5 Maximum number of iterations
                       as.integer(1),         #6 number of iteration
                       as.integer(1),         #7 Approximated depth
                       PACKAGE = "mrfDepth")

    deepest.fit <- result[[4]]
    deepest.fit <-deepest.fit[c(p1,1:(p1 - 1))]
    names(deepest.fit) <- c("intercept", paste("slope var.", 1:(p1 - 1)))

    returned.result <- list(deepest = deepest.fit,
                            depth = result[[7]] / n1,
                            ndir = result[[6]],
                            dimension = NULL,
                            hyperplane = NULL
                            )
    class(returned.result) <- c("mrfDepth", "rdepthmedian")
    return(returned.result)
  }
}
