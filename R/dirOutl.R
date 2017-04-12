

dirOutl <- function(x, z = NULL, options = list()){

  
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
    if(type == "compWise") {
      type <- 0
    }
    if(type == "Affine") {
      type <- 1
    }
    if(type == "Rotation") {
      type <- 2
    }
    if(type == "Shift") {
      type <- 3
    }

  } else {
    type <- 1 # Default "Affine"
  }
  if ("ndir"  %in% names(options)) {
    ndir <- options[["ndir"]]
  } else {
    if ( type == 1){
      ndir <- 250 * p1
    }else if (type == 2) {
      ndir <- 5000
    }else {
      ndir <- 12500
    }

  }
  if ("seed"  %in% names(options)) {
    seed <- options[["seed"]]
  } else {
    seed <- 10
  }
  if ("lb"  %in% names(options)) {
    lb   <- options[["lb"]]
  } else {
    lb   <- NULL
  }
  if ("precScale"  %in% names(options)) {
    precScale   <- options[["precScale"]]
  } else {
    precScale   <- 1e-10
  }
  if ("rmZeroes"  %in% names(options)) {
    rmZeroes   <- options[["rmZeroes"]]
  } else {
    rmZeroes   <- FALSE
  }
  if ("maxRatio"  %in% names(options)) {
    maxRatio   <- options[["maxRatio"]]
    if(maxRatio > 0 && maxRatio < 2){
      stop("maxRatio should be either 0 or larger than 2")
    }
  } else {
    maxRatio   <- 0
  }
  
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
    class(returned.result) <- c("mrfDepth", "dirOutl")
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
    class(returned.result) <- c("mrfDepth", "dirOutl")
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
    class(returned.result) <- c("mrfDepth", "dirOutl")
    return(returned.result)
  }
  set.seed(seed)
  res <- tryCatch( .Call("dirOutl_cpp", x, z, type, ndir, rmZeroes, maxRatio, precScale, 
                         PACKAGE = "mrfDepth"),
                   "std::range_error" = function(e){
                     conditionMessage( e ) })
  
  if( res$error_code == 1){
    warning("A direction was found for which the robust scale
              estimate equals\n zero. See the help page for more details.",
            call. = FALSE)
    
    hyperplane <- res$hyperplane
    datat <- x %*% hyperplane
    inSubspace <-  (abs(datat-median(datat)) < 1e-10)
    returned.result <- list(outlyingnessX = NULL, outlyingnessZ = NULL,
                            cutoff = NULL, flagX = NULL, flagZ = NULL,
                            singularSubsets = NULL, dimension = NULL,
                            hyperplane = hyperplane, inSubspace =inSubspace)
    class(returned.result) <- c("mrfDepth", "dirOutl")
    return(returned.result) 
  }
  
  outlyingnessX <- res$outlyingnessX
  outlyingnessZ <- res$outlyingnessZ
  
  cutoff <- sqrt(qchisq(0.99, p1)) * median(outlyingnessX)
  flag.X <- outlyingnessX <= cutoff
  flag.Z <- outlyingnessZ <= cutoff
  
  
  returned.result <- list(outlyingnessX = outlyingnessX,
                          outlyingnessZ = outlyingnessZ,
                          cutoff = cutoff,
                          flagX = flag.X,
                          flagZ = flag.Z,
                          singularsubsets = res$singularsubsets, 
                          dimension = NULL,
                          hyperplane = NULL,
                          inSubspace = NULL)
  class(returned.result) <- c("mrfDepth", "dirOutl")
  return(returned.result)
}