adjOutl <- function(x, z = NULL, options = list()){

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
  p1 <- ncol(x)
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

  #check options and load defaults if necessary.
  if (is.null(options)) {
    options <- list()
  }
  if (!is.list(options)) {
    stop("options must be a list")
  }
  if ("type" %in% names(options)) {
    type <- options[["type"]]
  } else {
    type <- "Affine"
  }
  if ("ndir"  %in% names(options)) {
    ndir <- options[["ndir"]]
  } else {
    ndir <- NULL
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
    if (type.id == 1) {
      ndir <- 250 * p1
    } else if (type.id == 2) {
      ndir <- 250 * 20
    } else {
      ndir <- 250 * 50
    }
  }
  calc.all <- 0
  # If the specified number of directions is larger than the
  # possible different directions, switch to exact computation.
  if (is.numeric(ndir)) {
    if (ndir < 1) {
      stop("The number of directions must be a positive integer.")
    }
    if (type.id == 1) {
      ndir0 <- choose(n1, p1)
      if (ndir0 <= ndir) {
        ndir <- ndir0
        calc.all <- 1
      }
    }
    if (type.id == 2) {
      ndir0 <- choose(n1, 2)
      if (ndir0 <= ndir) {
        ndir <- ndir0
        calc.all <- 1
      }
    }
  }
  if (!is.numeric(ndir)) {
    if (ndir == "all") {
      if (type.id == 1) {
        ndir <- choose(n1, p1)
        calc.all <- 1
        if (ndir > 1e7) {
          stop("ndir is larger than 1e7. Try a smaller value of ndir.")
        }
      }
      if (type.id == 2) {
        ndir <- choose(n1, 2)
        calc.all <- 1
        if (ndir > 1e7) {
          stop("ndir is larger than 1e7. Try a smaller value of ndir.")
        }
      }
      if (type.id == 3) {
        stop("Cannot compute all directions for type Shift.")
      }
    }
    else stop("The input parameter ndir is not recognized.")
  }
  if (type.id == 1 & (n1 < (p1 + 1))) {
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
    warning("One of the variables has zero
            standard deviation. Check the data matrix x.")
    returned.result <- list(outlyingnessX = NULL,
                            outlyingnessZ = NULL,
                            cutoff = NULL,
                            flagX = NULL,
                            flagZ = NULL,
                            singularSubsets = NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14),
                            inSubspace = NULL)
    class(returned.result) <- c("mrfDepth", "adjOutl")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n1 - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the
            data matrix x could not be computed.")
    returned.result <- list(outlyingnessX = NULL,
                            outlyingnessZ = NULL,
                            cutoff = NULL,
                            flagX = NULL,
                            flagZ = NULL,
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = NULL,
                            inSubspace = NULL)
    class(returned.result) <- c("mrfDepth", "adjOutl")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit is found. Check the output for more details.")
    returned.result <- list(outlyingnessX = NULL,
                            outlyingnessZ = NULL,
                            cutoff = NULL,
                            flagX = NULL,
                            flagZ = NULL,
                            singularSubsets = NULL,
                            dimension = sum(w1$d > tol),
                            hyperplane = w1$v[, which(w1$d == min(w1$d))[1]],
                            inSubspace = NULL)
    class(returned.result) <- c("mrfDepth", "adjOutl")
    return(returned.result)
  }

  #####
  # Perform the actual computations
  x <- rbind(x, z)
  n <- nrow(x)
  p <- ncol(x)

  result <- .C("adjprojout",
             as.integer(n),         #1 Total number of points.
             as.integer(p),         #2 Dimension of the data.
             as.integer(ndir),      #3 Number of directions.
             as.double(x),          #4 Data matrix (both x and z).
             as.double(rep(0, n)),  #5 Computed adjusted outlyingness.
             as.double(0),          #6 Medcouple computed on the adjusted
                                    #  outlying values of x.
             as.integer(0),         #7 Number of singular directions.
             as.integer(type.id),   #8 Integer indicating which type of
                                    #  directions to consider.
             as.integer(n1),        #9 Number of observations in x.
             as.integer(calc.all),  #10 Flag indicating whether all possible
                                    #   directions should be considered.
             as.double(rep(0, p1)), #11 Vector containing the direction on
                                    #   which a zero IQR is found.
             as.integer(seed),      #12 The seed.
             PACKAGE = "mrfDepth")

  adj.outlyingness <- result[[5]]
  if (sum(abs(result[[11]])) > tol) {
    warning("A direction was found for which the robust scale estimate equals
            zero. See the help page for more details.", call. = FALSE)
    returned.result <- list(outlyingnessX = NULL,
                            outlyingnessZ = NULL,
                            cutoff = NULL,
                            flagX = NULL,
                            flagZ = NULL,
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = result[[11]],
                            inSubspace = as.logical(adj.outlyingness[1:n1]))
    class(returned.result) <- c("mrfDepth", "adjOutl")
    return(returned.result)
  }

  LAO <- log(0.1 + adj.outlyingness[1:n1])
  cutoff <- exp(median(LAO) + mad(LAO) * qnorm(0.995)) - 0.1
  
  
  flag.X <- adj.outlyingness[1:n1] <= cutoff
  flag.Z <- adj.outlyingness[(n1 + 1):(n1 + n2)] <= cutoff

  returned.result <- list(outlyingnessX = adj.outlyingness[1:n1],
                          outlyingnessZ = adj.outlyingness[(n1 + 1):(n1 + n2)],
                          cutoff = cutoff,
                          flagX = flag.X,
                          flagZ = flag.Z,
                          singularsubsets = result[[7]],
                          dimension = NULL,
                          hyperplane = NULL,
                          inSubspace = NULL)
  class(returned.result) <- c("mrfDepth", "adjOutl")
  return(returned.result)
}
