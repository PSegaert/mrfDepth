projdepth <- function(x, z = NULL, options = list()) {

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
    stop("The data dimension has to be the same for x and z.")
  }
  if (n2 > sum(complete.cases(z))) {
    stop("Missing values in z are not allowed.")
  }
  #check options
  if (is.null(options)) {
    options <- list()
  }
  if (!is.list(options)) {
    stop("options must be a list")
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
    returned.result <- list(depthX = NULL,
                            depthZ = NULL,
                            cutoff = NULL,
                            flagX = NULL,
                            flagY = NULL,
                            singularSubsets  =  NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14),
                            inSubspace = NULL)
    class(returned.result) <- c("mrfDepth", "projdepth")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n1 - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix
            x could not be computed.")
    returned.result <- list(depthX = NULL,
                            depthZ = NULL,
                            cutoff = NULL,
                            flagX = NULL,
                            flagY = NULL,
                            singularSubsets  =  NULL,
                            dimension = NULL,
                            hyperplane = NULL,
                            inSubspace = NULL)
    class(returned.result) <- c("mrfDepth", "projdepth")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more details.")
    returned.result <- list(depthX = NULL,
                            depthZ = NULL,
                            cutoff = NULL,
                            flagX  = NULL,
                            flagZ = NULL,
                            singularSubsets = NULL,
                            dimension = sum(w1$d > tol),
                            hyperplane = w1$v[which(w1$d == min(w1$d))[1]],
                            inSubspace = NULL)
    class(returned.result) <- c("mrfDepth", "projdepth")
    return(returned.result)
  }

  original <- options(warn = 1)
  result <- outlyingness(x = x, z = z, options = options)
  options(warn = original$warn)

  if (!is.null(result$hyperplane)) {
    returned.result <- list(depthX = NULL,
                            depthZ = NULL,
                            cutoff = NULL,
                            flagX = NULL,
                            flagZ = NULL,
                            singularSubsets = NULL,
                            dimension = NULL,
                            hyperplane = result[["hyperplane"]],
                            inSubspace = result[["inSubspace"]])
    class(returned.result) <- c("mrfDepth", "projdepth")
    return(returned.result)
  }
  else{
    returned.result <- list(depthX = 1 / (1 + result[["outlyingnessX"]]),
                            depthZ = 1 / (1 + result[["outlyingnessZ"]]),
                            cutoff = 1 / (1 + result[["cutoff"]]),
                            flagX = result[["flagX"]],
                            flagZ = result[["flagZ"]],
                            singularSubsets = result[["singularSubsets"]],
                            dimension = NULL,
                            hyperplane = result[["hyperplane"]],
                            inSubspace = NULL)
    class(returned.result) <- c("mrfDepth", "projdepth")
    return(returned.result)
  }

}
