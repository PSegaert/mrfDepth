fOutl <- function(x,
                  z = NULL,
                  type = "fAO",
                  alpha = 0,
                  time = NULL,
                  diagnostic = FALSE,
                  distOptions=NULL) {

  ######
  # Check input.
  if (missing(x)) {
    stop("Input argument x is required.")
  }

  #Check x
  if (!is.array(x)) {
    stop("x must be a three dimensional array.")
  }
  if (length(dim(x)) != 3) {
    stop("x must be a three dimensional array.")
  }
  if (sum(is.nan(x)) != 0) {
    stop("x contains missing cases.")
  }
  t1 <- dim(x)[1]
  n1 <- dim(x)[2]
  p1 <- dim(x)[3]

  #Check z
  if (is.null(z)) {
    z <- x
  }
  if (!is.array(z)) {
    stop("z must be a three dimensional array.")
  }
  if (length(dim(z)) != 3) {
    stop("z must be a three dimensional array.")
  }
  if (sum(is.nan(z)) != 0) {
    stop("z contains missing cases.")
  }
  t2 <- dim(z)[1]
  n2 <- dim(z)[2]
  p2 <- dim(z)[3]

  #Check dimension match between x and z
  if (p1 != p2) {
    stop("The p dimension of x and z must match.")
  }
  if (t1 != t2) {
    stop("The t dimension of x and z must match.")
  }

  #Check type
  Indtype <- match(type, c("fSDO", "fAO",
                           "fDO", "fbd"))[1]
  if (is.na(Indtype)) {
    stop("type should be one of fSDO, fAO, fDO or fbd.")
  }

  #Check alpha
  if (!is.numeric(alpha)) {
    stop("alpha must be numeric")
  }
  else if (is.vector(alpha)) {
   if(alpha != 0) {
     stop("Only the value 0 for alpha is allowed.")
   }
  }
  else if (is.matrix(alpha)) {
    NRowAlpha <- dim(alpha)[1]
    NColAlpha <- dim(alpha)[2]
    if (NRowAlpha != 1 || NColAlpha != t1) {
      stop("alpha must be a (1xt)-row matrix.")
    }
  }
  else{
    stop("alpha must be either a number or a (1xt)-row matrix.")
  }

  #Check time
  if (is.null(time)) {
    time <- 1:t1
  }
  if (!is.numeric(time) || !is.vector(time)) {
    stop("time should be a numeric vector.")
  }
  if (length(time) != t1) {
    stop("time should contain t elements")
  }
  if (length(time) != 1) {
    dTime <- diff(c(time[1], time, time[t1]), lag = 2)
    if (min(dTime) <= 0) {
      stop("time should be strictly increasing.")
    }
  }
  else{
    dTime <- 1
  }

  #check diagnostic
  if (!is.logical(diagnostic)) {
    stop("diagnostic should be a logical")
  }

  #check depthOptions
  if (is.null(distOptions)) {
    distOptions <- list()
  }
  if (!is.list(distOptions)) {
    stop("distOptions must be a list")
  }


  weights <- rep(1, t1)
  distTimeX <- matrix(NA, nrow = n1, ncol = t1)
  distTimeZ <- matrix(0.0, nrow = n2, ncol = t2)
  locOutlX <- matrix(NA, nrow = n1, ncol = t1)
  locOutlZ <- matrix(0, nrow = n2, ncol = t2)
  if (is.matrix(alpha)) {
    weights <- alpha
  }

  warningFlagFit <- 0
  warningIndFit <- c()

  Original <- options(warn = -1)
  for (j in 1:t1) {
    exactfit <- 0

    #R has standard dimension dropping, we need to be carefull
    if (p1 == 1)  {
      xTimePoint <- matrix(x[j,,1])
      zTimePoint <- matrix(z[j,,1])
    }
    else{
      xTimePoint <- x[j,,,drop = TRUE]
      zTimePoint <- z[j,,,drop = TRUE]
    }

    #Find cross-sectional distance
    if (type == "fSDO") {
      temp <- outlyingness(x = xTimePoint,
                           z = zTimePoint,
                           options  =  distOptions)
      if (!is.list(temp)) {
        temp <- list()
      }
      if (!is.null(temp$outlyingnessZ)) {
        distTimeX[,j] <- temp$outlyingnessX
        distTimeZ[,j] <- temp$outlyingnessZ
        locOutlX[,j] <- as.numeric(!temp$flagX)
        locOutlZ[,j] <- as.numeric(!temp$flagZ)
      }
      else{
        exactfit <- 1
      }
    }
    else if (type == "fAO") {
      temp <- adjOutl(x = xTimePoint,
                      z = zTimePoint,
                      options = distOptions)
      if (!is.null(temp$outlyingnessZ)) {
        distTimeX[,j] <- temp$outlyingnessX
        distTimeZ[,j] <- temp$outlyingnessZ
        locOutlX[,j] <- as.numeric(!temp$flagX)
        locOutlZ[,j] <- as.numeric(!temp$flagZ)
      }
      else{
        exactfit <- 1
      }
    }
    else if (type == "fDO") {
      temp <- dirOutl(x = xTimePoint,
                      z = zTimePoint,
                      options = distOptions)
      if (!is.null(temp$outlyingnessZ)) {
        distTimeX[,j] <- temp$outlyingnessX
        distTimeZ[,j] <- temp$outlyingnessZ
        locOutlX[,j] <- as.numeric(!temp$flagX)
        locOutlZ[,j] <- as.numeric(!temp$flagZ)
      }
      else{
        exactfit <- 1
      }
    }
    else if (type == "fbd") {
      temp <- bagdistance(x = xTimePoint,
                          z = zTimePoint,
                          options = distOptions)
      if (!is.null(temp$bagdistance)) {
        distTimeZ[,j] <- temp$bagdistance
        locOutlZ[,j] <- as.numeric(!temp$flag)
      }
      else{
        exactfit <- 1
      }
    }

    #Check if exact fit needs handling later on
    if (exactfit) {
      weights[j] <- 0
      warningFlagFit <- 1
      warningIndFit <- c(warningIndFit, j)
    }

  }

  options(warn = Original$warn)

  weights <- weights * dTime
  weights <- weights / sum(weights)
  fOutlX <- distTimeX %*% weights
  fOutlZ <- distTimeZ %*% weights

  #Assemble the results
  Result <- list(fOutlyingnessX = fOutlX,
                 fOutlyingnessZ = fOutlZ,
                 weights = weights)
  if (diagnostic) {
    Result$crossDistsX <- distTimeX
    Result$crossDistsZ <- distTimeZ
    Result$locOutlX <- locOutlX
    Result$locOutlZ <- locOutlZ
  }
  Result$distType <- type
  class(Result) <- c("mrfDepth", "fOutl")

  #Handle all warnings
  if (warningFlagFit == 1) {
    warning(paste("Exact fits were detected for certain time points.",
                  "Their weights will be set to zero.",
                  "Check the returned results"))
    Result$IndFlagExactFit <- warningIndFit
  }

  return(Result)
}
