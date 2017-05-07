mfd <- function(x,
              z = NULL,
              type = "hdepth",
              alpha = 0,
              time = NULL,
              diagnostic = FALSE,
              depthOptions = NULL) {

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
  Indtype <- match(type, c("hdepth", "projdepth",
                           "sprojdepth", "dprojdepth", "sdepth"))[1]
  if (is.na(Indtype)) {
    stop("type should be one of hdepth, projdepth , sprojdepth, dprojdepth or sdepth.")
  }
  if (Indtype == 5 && p1 > 2) {
    stop("sdepth depth only implemented for p<=2.")
  }

  #Check alpha
  if (!is.numeric(alpha)) {
    stop("alpha must be numeric")
  }
  if (is.vector(alpha)) {
    if (alpha < 0) {
      stop("alpha should be part of [0,1]")
    }
    if (alpha > 1) {
      stop("alpha should be part of [0,1]")
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
  if (is.null(depthOptions)) {
    depthOptions <- list()
  }
  if (!is.list(depthOptions)) {
    stop("depthOptions must be a list")
  }


  weights <- rep(1, t1)
  depthsTimeX <- matrix(NA, nrow = n1, ncol = t1)
  depthsTimeZ <- matrix(0.0, nrow = n2, ncol = t2)
  locOutlX <- matrix(NA, nrow = n1, ncol = t1)
  locOutlZ <- matrix(NA, nrow = n1, ncol = t1)
  if (is.matrix(alpha)) {
    weights <- alpha
  }

  warningFlagFit <- warningFlagBag <- warningFlagIso <- warningFlagAlpha <- 0
  warningIndFit <- warningIndBag <- warningIndIso <- warningIndAlpha <- c()

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

    #Find cross-sectional depth
    if (type == "hdepth") {
      temp <- hdepth(x = xTimePoint, z = zTimePoint, options  =  depthOptions)
      if (!is.list(temp)) {
        temp <- list()
      }
      if (!is.null(temp$depthZ)) {
        depthsTimeX[,j] <- temp$depthX
        depthsTimeZ[,j] <- temp$depthZ
      }
      else{
        exactfit <- 1
      }

      #If requested find local halfspace outliers
      if (diagnostic & p1 == 2 & exactfit == FALSE) {
          temp <- compBagplot(x = xTimePoint, type = type)
          if (sum(is.nan(temp$flag)) == 0) {
            locOutlX[,j] <- temp$flag
          }
          else{
            warningFlagBag <- 1
            warningIndBag <- c(warningIndBag, j)
          }
      }

    }
    else if (type == "projdepth") {
      temp <- projdepth(x = xTimePoint, z = zTimePoint, options = depthOptions)
      if (!is.null(temp$depthZ)) {
        depthsTimeX[,j] <- temp$depthX
        depthsTimeZ[,j] <- temp$depthZ
        locOutlX[,j] <- as.numeric(!temp$flagX)
        locOutlZ[,j] <- as.numeric(!temp$flagZ)
      }
      else{
        exactfit <- 1
      }
    }
    else if (type == "sprojdepth") {
      temp <- sprojdepth(x = xTimePoint, z = zTimePoint, options = depthOptions)
      if (!is.null(temp$depthZ)) {
        depthsTimeX[,j] <- temp$depthX
        depthsTimeZ[,j] <- temp$depthZ
        locOutlX[,j] <- as.numeric(!temp$flagX)
        locOutlZ[,j] <- as.numeric(!temp$flagZ)
      }
      else{
        exactfit <- 1
      }
    }
    else if (type == "dprojdepth") {
      temp <- dprojdepth(x = xTimePoint, z = zTimePoint, options = depthOptions)
      if (!is.null(temp$depthZ)) {
        depthsTimeX[,j] <- temp$depthX
        depthsTimeZ[,j] <- temp$depthZ
        locOutlX[,j] <- as.numeric(!temp$flagX)
        locOutlZ[,j] <- as.numeric(!temp$flagZ)
      }
      else{
        exactfit <- 1
      }
    }
    else{
      temp <- sdepth(x = xTimePoint, z = zTimePoint)
      if (!is.null(temp$depth)) {
        depthsTimeZ[,j] <- temp$depth
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

    #Calculate the area of depth region at time T
    #Do not calculate when alpha is row-matrix
    if (!is.matrix(alpha)) {
      #Only for non-constant weights, no point in calculating for exact fits
      if (alpha != 0 && exactfit == 0) {
        temp <- depthContour(x = xTimePoint, alpha, type = type)
        Vert <- temp[[1]]$vertices
        if (sum(is.nan(Vert)) == 0) {
          if (nrow(Vert) != nrow(unique(Vert))) {
            warningFlagIso <- 1
            warningIndIso <- c(warningIndIso, j)
          }
          else{
            if (p1 == 1) {
              temp <- max(Vert) - min(Vert)
            } else {
              temp <- try(convhulln(matrix(Vert, ncol = p1), "FA")$vol,
                          silent = TRUE)
            }
            if (!is.numeric(temp)) {
              warningFlagAlpha <- 1
              warningIndAlpha <- c(warningIndAlpha, j)
            }
            else{
              weights[j] <- temp
            }
          }
        }
        else{
          weights[j] <- 0
          warningFlagIso <- 1
          warningIndIso <- c(warningIndIso, j)
        }
      }
    }
  }

  options(warn = Original$warn)

  weights <- weights * dTime
  weights <- weights / sum(weights)
  depthsX <- depthsTimeX %*% weights
  depthsZ <- depthsTimeZ %*% weights

  #Assemble the results
  Result <- list(MFDdepthX = depthsX,
                 MFDdepthZ = depthsZ,
                 weights = weights)
  if (diagnostic) {
    Result$crossdepthX <- depthsTimeX
    Result$crossdepthZ <- depthsTimeZ
    Result$locOutlX <- locOutlX
    Result$locOutlZ <- locOutlZ
  }
  Result$depthType <- type
  class(Result) <- c("mrfDepth", "mfd")

  #Handle all warnings
  if (warningFlagFit == 1) {
    warning(paste("Exact fits were detected for certain time points.",
                  "Their weights will be set to zero.",
                  "Check the returned results"),
            call. = FALSE)
    Result$IndFlagExactFit <- warningIndFit
  }
  if (warningFlagBag == 1) {
    warning(paste("The bagplot could not be computed at all time points.",
                  "Their weights will be set to zero.",
                  "Check the returned results"),
            call. = FALSE)
    Result$IndFlagBag <- warningIndBag
  }
  if (warningFlagIso == 1) {
    warning(paste("The isohdepth contours could not be computed at all", 
                  "time points. Their weights will be set to zero.",
                  "Check the returned results"),
            call. = FALSE)
    Result$IndFlagIso <- warningIndIso
  }
  if (warningFlagAlpha  ==  1) {
    warning(paste("The specified alpha is too large at all time points.",
                  "Their weights will be set to zero.",
                  "Check the returned results"),
            call. = FALSE)
    Result$IndFlagAlpha <- warningIndAlpha
  }

  return(Result)
}
