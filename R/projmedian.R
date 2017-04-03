projmedian <- function(x, projection.depths = NULL, options = NULL){

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
  #check projection.depths
  if (!is.null(projection.depths)) {
    projection.depths <- data.matrix(projection.depths)
    n2 <- nrow(projection.depths)
    p2 <- ncol(projection.depths)
    if (n1 != n2) {
      stop("A different number of depths from the number of
           observations was specified.")
    }
    if (p2 != 1) {
      stop("Provided depths have to be a columnvector.")
    }
    if (sum(projection.depths > 1) != 0) {
      stop("The user supplied depths must take values in ]0,1].")
    }
    if (sum(projection.depths <= 0) != 0) {
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


  if (is.null(projection.depths)) {
    return.depth.calculation <- projdepth(x = x, options = options)
    projection.depths <- return.depth.calculation[["depthX"]]
  }

  returned.result <- findCenterProj(x, projection.depths)
  class(returned.result) <- c("mrfDepth", "projmedian")
  
  return(returned.result)

}


findCenterProj <- function(x, projection.depths){
  n <- nrow(x)
  p <- ncol(x)
  center <- vector("list", 4)

  #max.depth
  max.depth <- max(projection.depths)
  ind.max.depths <- which(projection.depths == max.depth)
  center[[1]] <- colMeans(matrix(x[ind.max.depths, ], ncol = p))
  center[[1]] <- matrix(center[[1]], ncol = p)
  names(center)[1] <- "max"


  #Gravity
  cut.off <- median(projection.depths)
  ind.grav.depths <- which(projection.depths >= cut.off)
  center[[2]] <- colMeans(matrix(x[ind.grav.depths, ], ncol = p))
  center[[2]] <- matrix(center[[2]], ncol = p)
  names(center)[2] <- "gravity"


  #Huber
  weights <- matrix(rep(1, n), ncol = 1)
  cut.off <- sqrt(qchisq(0.95, df = p))
  ind.huber <- which( projection.depths <= (1 / (1 + cut.off)) )
  outlyingness <- 1 / projection.depths[ind.huber] - 1
  weights[ind.huber] <- (cut.off / outlyingness) ^ 2
  center[[3]] <- t( (t(x) %*% weights) / n )
  center[[3]] <- matrix(center[[3]], ncol = p)
  names(center)[3] <- "Huber"

  return(center)
}
