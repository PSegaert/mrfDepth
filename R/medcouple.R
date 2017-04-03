medcouple <- function(x, do.reflect=NULL){

  ######
  # Check input.
  if (missing(x)) {
    stop("Input argument x is required.")
  }

  #Check the x data.
  x <- data.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if (n > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }
  if (!is.numeric(x)) {
    stop("x should be a numeric data matrix.")
  }
  #check reflection
  if (is.null(do.reflect)) {
    if (n > 100) {
      do.reflect <- TRUE
    }
    else {
      do.reflect <- FALSE
    }
  }
  else{
    if (!(do.reflect %in% c(FALSE, TRUE))) {
      stop("doreflect should be one of TRUE or FALSE.")
    }
  }

  mc.result <- rep(0.0, p)

  for (i in 1:p) {

    temp <- .C("medcoupleC",
              as.double(x[, i]),     #1 Data vector
              as.integer(n),         #2 Number of observations
              as.double(0.0),        #3 Medcouple
              as.integer(do.reflect), #4 Logical indicating calculation on -x
              PACKAGE = "mrfDepth")

    mc.result[i] <- temp[[3]]

  }

  class(mc.result) <- c("mrfDepth", "medcouple")
  
  return(mc.result)

}
