sprojmedian <- function(x, sprojection.depths = NULL, options = NULL) {

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
  if (n1 > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }
  #check spDepths
  if (!is.null(sprojection.depths)) {
    sprojection.depths <- data.matrix(sprojection.depths)
    n2 <- nrow(sprojection.depths)
    p2 <- ncol(sprojection.depths)
    if (n1 != n2) {
      stop("A different number of depths from the number of
           observations was specified.")
    }
    if (p2 != 1) {
      stop("Provided depths have to be a columnvector.")
    }
    if (sum(sprojection.depths > 1) != 0) {
      stop("The user supplied depths must take values in ]0,1].")
    }
    if (sum(sprojection.depths <= 0) != 0) {
      stop("The user supplied depths must take values in ]0,1].")
    }
  }
  #check options
  if (is.null(options)) {
    options <- list()
  }
  if (!is.list(options)) {
    stop("options must be a list")
  }

  if (is.null(sprojection.depths)) {
    return.depth.calculation <- sprojdepth(x = x, options = options)
    sprojection.depths <- return.depth.calculation[["depthX"]]
  }

  returned.result <- findCenterSProj(x, sprojection.depths)
  class(returned.result) <- c("mrfDepth", "sprojmedian")
  
  return(returned.result)

}


findCenterSProj <- function(x, sprojection.depths) {
  p <- ncol(x)
  center <- vector("list", 2)

  #MaxDepth
  max.depth <- max(sprojection.depths)
  ind.max.depths <- which(sprojection.depths == max.depth)
  center[[1]] <- colMeans(matrix(x[ind.max.depths, ], ncol = p))
  center[[1]] <- matrix(center[[1]], ncol = p)
  names(center)[1] <- "max"

  #Gravity
  cut.off <- median(sprojection.depths)
  ind.grav.depths <- which(sprojection.depths >= cut.off)
  center[[2]] <- colMeans(matrix(x[ind.grav.depths, ], ncol = p))
  center[[2]] <- matrix(center[[2]], ncol = p)
  names(center)[2] <- "gravity"

  return(center)
}
