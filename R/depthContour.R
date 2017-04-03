depthContour <- function(x, alpha = NULL, type = "hdepth", directions = NULL,
                         options=NULL){

  ######
  # Check input.
  if (missing(x)) {
    stop("The input argument x is required.")
  }

  # Check the x data.
  x <- data.matrix(x)
  if (sum(is.nan(x) != 0)) {
    stop("Missing values in matrix x are not allowed.")
  }
  if (!is.numeric(x)) {
    stop("The input argument x must be a numeric data matrix.")
  }
  n <- nrow(x)
  p <- ncol(x)
  if (n < (p + 1)) {
    stop("The number of observations should be larger than the number
         of variables + 1.")
  }

  # check alpha
  if (!is.numeric(alpha)) {
    stop("alpha should be a numeric vector.")
   }
  if (sum(sort(alpha) == alpha) != length(alpha)) {
    stop("The vector alpha should be sorted in ascending order.")
  }

  # Check type
  Indtype <- match(type, c("hdepth", "projdepth", "sprojdepth"))[1]
  if (is.na(Indtype)) {
    stop("type should be one of hdepth, projdepth or sprojdepth.")
  }
  if (min(alpha) < 1 / n && Indtype == 1) {
    stop("alpha should be at least 1/n for halfspace depth.")
  }
  if (max(alpha) > 0.5 && Indtype == 1) {
    stop("alpha should be smaller than 0.5.")
  }
  if (min(alpha) <= 0 || max(alpha) >= 1) {
    stop("alpha should lie in the open interval (0,1).")
  }

  # Check directions
  if (is.null(directions)) {
      set.seed(123)
      NDirections <- 250 * p
      directions <- matrix(rnorm(n = NDirections * p), ncol = p)
      if (p == 1) {
        directions <- matrix(c(-1,1), ncol = 1)
      }
      directions <- directions / sqrt(rowSums(directions ^ 2))
  } else{
    directions <- data.matrix(directions)
    if (sum(is.nan(directions) != 0)) {
      stop("Missing values in the matrix directions are not allowed.")
    }
    if (!is.numeric(directions)) {
      stop("The input argument directions must be a numeric data matrix.")
    }
    p.dir <- ncol(directions)
    if (p.dir != p) {
      stop("The dimension of the matrix directions should be the same as
           the dimension of x.")
    }
  }
  directions <- directions / sqrt(rowSums(directions ^ 2))

  # check options
  if (is.null(options)) {
    options <- list(approx = TRUE)
  }
  if (!is.list(options)) {
    stop("The input argument options must be a list")
  }
  if ("max.iter"  %in% names(options)) {
    max.iter <- options[["max.iter"]]
    if (is.numeric(max.iter)) {
      if (max.iter < 1) {
        stop("The option max.iter must be a strictly positive integer.")
      }
    }
    if (!is.numeric(max.iter)) {
      stop("The option max.iter must be a strictly positive integer.")
    }
  } else {
    options$max.iter <- 100
  }


  # Check data for exact fits
  tol <- 1e-7
  scaled.x <- scale(x)
  temp <- attributes(scaled.x)
  column.sd <- temp[["scaled:scale"]]
  if (sum(column.sd <= 1e-14) > 0) {
    warning("At least one of the variables has zero standard deviation.
            Check the data matrix x.")
    Result <- list(depth = NULL,
                   vertices = NULL,
                   empty = NULL,
                   dithered = NULL,
                   converged = NULL,
                   type = type,
                   dimension = sum(column.sd > 1e-14),
                   hyperplane = as.numeric(column.sd <= 1e-14))
    class(Result) <- c("mrfDepth", "depthContour")
    return(Result)
  }
  w1 <- try(svd(scaled.x / sqrt(n - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x could
            not be computed.")
    Result <- list(depth = NULL,
                 vertices = NULL,
                 empty = NULL,
                 dithered = NULL,
                 converged = NULL,
                 type = type,
                 dimensions = NULL,
                 hyperplane = NULL)
    class(Result) <- c("mrfDepth", "depthContour")
    return(Result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more details.")
    Result <- list(depth = NULL,
                   vertices = NULL,
                   empty = NULL,
                   dithered = NULL,
                   converged = NULL,
                   type = type,
                   dimensions = sum(w1$d > tol),
                   hyperplane = w1$v[, which(w1$d == min(w1$d))[1]])
    class(Result) <- c("mrfDepth", "depthContour")
    return(Result)
  }

  if (type == "hdepth") {
    if (p == 2) {
      Result <- CalcBivHContour(x, alpha)
    }
    if (p != 2) {
      Result <- CalcHContour(x, alpha, directions, options)
    }
  }
  else if (type == "projdepth" || type == "sprojdepth") {
    Result <- CalcProjTypeContour(x, alpha, type, directions, options)
  }
  else{
    stop("The input argument type is not recognized.")
  }

  class(Result) <- c("mrfDepth", "depthContour")
  return(Result)
}

CalcBivHContour <- function(x, alpha) {

  n <- nrow(x)

  Result <- vector("list", length(alpha) + 1)
  names(Result)[[length(alpha) + 1]] <- "type"

  for (i in 1:length(alpha)) {
    TResult <- .Fortran("iso2hdw",
      as.double(x[, 1, drop = TRUE]),  #1 First variable of the data set.
      as.double(x[, 2, drop = TRUE]),  #2 Second variable of the data set.
      as.integer(n),                   #3 Number of observations in the data set
      as.integer(floor(alpha[i] * n)), #4 Depth of the contour to be calculated.
      as.double(rep(0, n * (n - 1) / 2)),  #5 First coordinate of the vertices.
      as.double(rep(0, n * (n - 1) / 2)),  #6 Second coordinate of the vertices.
      as.integer(1),                   #7 Logical signaling empty contour
      as.integer(1),                   #8 Number of vertices in the contour.
      as.integer(1),                   #9 Logical signaling whether dithering
                                       #  was performed.
      PACKAGE = "mrfDepth")

    names(Result)[[i]] <- "Contour"
    Result[[i]] <- list(depth = TResult[[4]] / n,
                        vertices = cbind(TResult[[5]][1:TResult[[8]]],
                                         TResult[[6]][1:TResult[[8]]]),
                        empty = as.logical(TResult[[7]]),
                        dithered = as.logical(TResult[[9]]),
                        converged = rep(TRUE, n),
                        type = "hdepth",
                        dimensions = NULL,
                        hyperplane = NULL)
  }

  Result[[length(alpha) + 1]] <- "hdepth"
  class(Result) <- c("mrfDepth", "depthContour")
  return(Result)

}

CalcHContour <- function(x, alpha, directions, options) {

  NObs <- nrow(x)

  TRes <- hdepthmedian(x = x)
  Center <- TRes$median
  center.depth <- TRes$depth

  Result <- vector("list", length(alpha) + 1)
  names(Result)[[length(alpha) + 1]] <- "type"

  for (i in 1:length(alpha)) {
    if (floor(NObs * alpha[i]) / NObs > center.depth) {
      Result[[i]] <- list(depth = floor(NObs * alpha[i]) / NObs,
                          vertices = NULL,
                          empty = TRUE,
                          dithered = FALSE,
                          converged = NULL)
    } else {
      cont.result <- CalcOneHalfContourIntersect(x, Center, alpha[i],
                                                directions, options)
      names(Result)[[i]] <- "Contour"
      Result[[i]] <- list(depth = floor(NObs * alpha[i]) / NObs,
                          vertices = cont.result$vertices,
                          empty = FALSE,
                          dithered = FALSE,
                          converged = cont.result$converged,
                          type = "hdepth",
                          dimensions = NULL,
                          hyperplane = NULL)
    }
  }

  Result[[length(alpha) + 1]] <- "hdepth"
  class(Result) <- c("mrfDepth", "depthContour")
  return(Result)

}

CalcOneHalfContourIntersect <- function(x, Center, alpha, Directions, options) {

  NObs <- nrow(x)
  NVar <- ncol(x)
  NPoints <- nrow(Directions)

  #Initialise the depth to target
  target.depth <- floor(NObs * alpha) / NObs

  #Center the data
  distribution.data <- sweep(x, MARGIN = 2, Center, FUN = "-")

  #Find the upper bounds by projection onto the direction.
  #The multivariate depth of points on the contour is smaller
  #than the univariate depth on the projection. Therefore
  #we can take the point with univariate depth equal
  #to the target.depth as an upper limit.
  Upperbounds <- matrix(0.0, nrow = NPoints, ncol = NVar)
  for (i in 1:NPoints) {
    projection.data <- distribution.data %*% Directions[i, ]
    #The vector in the matrix Distribution has a direction defined.
    #Therefore the target point must be situated in the second half
    #of the projection.
    temp.ind <- NObs - ceiling(target.depth * NObs)
    start.dev <- sort(projection.data, partial = temp.ind)[temp.ind]
    Upperbounds[i, ] <- 1.2 * abs(start.dev) * Directions[i, ]
  }

  #Lets go
  Lowerbounds <- matrix(0.0, ncol = NVar, nrow = NPoints)
  considered.points <- 0.5 * Upperbounds + Lowerbounds

  #Initialise constant for iterative algorithm
  maxit <- options$max.iter
  count <- 1
  ind.improve <- rep(TRUE, NPoints)
  d.cons.points <- rep(0.0, NPoints)
  Converged <- rep(TRUE, NPoints)

  #Exception for points equal to the median
  ind.cancel <- which(rowSums(Upperbounds) == 0)
  ind.improve[ind.cancel] <- FALSE

  Directions <- NULL
  while ((sum(ind.improve == TRUE) != 0) & (count <= maxit) ) {
    Result <- hdepth(x = distribution.data,
                     z = considered.points[ind.improve, , drop = FALSE],
                     options = options)
    d.cons.points[ind.improve] <- Result[["depthZ"]]
    if (is.null(Directions)) {
      Directions <- Result$direction
    }
    #Find out which are close enough
    ind.stop <- which(abs(d.cons.points - target.depth) <= 1 / (2 * NObs))
    ind.improve[ind.stop] <- FALSE
    d.cons.points[ind.stop] <- target.depth
    #Find out which need to be closer to the center
    ind.closer <- which(d.cons.points < target.depth)
    Upperbounds[ind.closer, ] <- considered.points[ind.closer, ]
    considered.points[ind.closer, ] <- 0.5 * (Lowerbounds[ind.closer, ] +
                                              considered.points[ind.closer, ])
    #Find out which need to be further away from the center
    ind.further <- which(d.cons.points > target.depth)
    Lowerbounds[ind.further, ] <- considered.points[ind.further, ]
    considered.points[ind.further, ] <- 0.5 * (considered.points[ind.further,
                                                                 ] +
                                               Upperbounds[ind.further, ])

    #Add a breaking condition for extreme slow convergence situations
    temp1 <- sqrt(rowSums(matrix(abs(Lowerbounds - Upperbounds),
                                 ncol = NVar) ^ 2
                          )
                  )
    temp2 <- sqrt(rowSums(Upperbounds ^ 2))
    aid.ind <- which(temp1 / temp2 <= 10 ^ -10)
    ind.improve[aid.ind] <- FALSE

    count <- count + 1
  }

  considered.points <- sweep(x = considered.points,
                             MARGIN = 2, Center, FUN = "+")
  Converged[ind.improve] <- FALSE

  return(list(vertices = considered.points,
              converged = Converged)
  )
}

CalcProjTypeContour <- function(x, alpha, type, Directions, options) {

  if (type == "projdepth") {
    TResult <- projdepth(x = x, z = x, options = options)
    Center <- projmedian(x = x,
                         projection.depths = TResult[["depthZ"]],
                         options = options)$max
  }
  else if (type == "sprojdepth") {
    TResult <- sprojdepth(x = x, z = x, options = options)
    Center <- sprojmedian(x = x,
                          sprojection.depths = TResult[["depthZ"]],
                          options = options)$max
  }
  else{
    stop("The input parameter type is not recognized.")
  }

  Result <- vector("list", length(alpha) + 1)
  names(Result)[[length(alpha) + 1]] <- "type"

  for (i in 1:length(alpha)) {
    if (alpha[i] > max(TResult[["depthX"]])) {
      Result[[i]] <- list(depth = alpha[i],
                          vertices = NULL,
                          empty = TRUE,
                          dithered = FALSE,
                          converged = NULL,
                          type = type,
                          dimensions = NULL,
                          hyperplane = NULL)
    } else {
      cont.result <- CalculateOneProjTypeContour(x, Center, alpha[i], type,
                                                 Directions, options)
      names(Result)[[i]] <- "Contour"
      Result[[i]] <- list(depth = alpha[i],
                          vertices = cont.result$vertices,
                          empty = FALSE,
                          dithered = FALSE,
                          converged = cont.result$converged,
                          type = type,
                          dimensions = NULL,
                          hyperplane = NULL)
    }
  }

  Result[[length(alpha) + 1]] <- type
  class(Result) <- c("mrfDepth", "depthContour")
  return(Result)

}

CalculateOneProjTypeContour <- function(x, Center, alpha, type,
                                        Directions, options) {

  NDirections <- nrow(Directions)
  NVar <- ncol(x)

  #Center the data
  x <- sweep(x, MARGIN = 2, Center, FUN = "-")

  #Initialise target depth
  target.depth <- alpha
  target.outlying <- (1 / target.depth) - 1

  #Find the upper bounds
  Upperbounds <- matrix(0.0, nrow = NDirections, ncol = NVar, byrow = TRUE)
  for (i in 1:NDirections) {
    projection.data <- x %*% Directions[i, ]
    if (type == "projdepth") {
      Scale <- mad(projection.data)
      Location <- median(projection.data)
      Upperbounds[i, ] <- (target.outlying * Scale + Location) * Directions[i, ]
    }
    else if (type == "sprojdepth") {
      Quants <- quantile(projection.data, probs = c(0.25, 0.5,  0.75), type = 5)
      IQR <- Quants[3] - Quants[1]
      Location <- Quants[2]
      MC <- medcouple(projection.data)
      if (MC > 0) {
        MC <- 3 * MC
      } else {
        MC <- 4 * MC
      }
      w2 <- Quants[3] + 1.5 * exp(MC) * IQR
      Upperbounds[i, ] <- 1.2 * (target.outlying *
                                  (w2 - Location) + Location ) * Directions[i, ]
    }
    else{
      stop("Type not recognized")
    }
  }

  Lowerbounds <- matrix(0.0, nrow = NDirections, ncol = NVar, byrow = TRUE)

  # Initialise starting points
  considered.points <- 0.5 * (Upperbounds + Lowerbounds)

  maxit <- options$max.iter
  count <- 1
  ind.improve <- rep(TRUE, NDirections)
  d.cons.points <- rep(0.0, NDirections)
  Converged <- rep(TRUE, NDirections)

  while ((sum(ind.improve == TRUE) != 0) & (count <= maxit) ) {
    if (type == "projdepth") {
      d.cons.points[ind.improve] <- projdepth(x = x,
                               z = considered.points[ind.improve, ,
                                                     drop = FALSE],
                               options = options)$depthZ
    }
    else if (type == "sprojdepth") {
      d.cons.points[ind.improve] <- sprojdepth(x = x,
                              z = considered.points[ind.improve, ,
                                                    drop = FALSE],
                              options = options)$depthZ
    }
    else{
      stop("Type not recognized.")
    }

    #Find out which are close enough
    ind.stop <- which(abs(d.cons.points - target.depth) <=  10 ^ -3)
    ind.improve[ind.stop] <- FALSE
    d.cons.points[ind.stop] <- target.depth
    #Find out which need to be closer to the center
    ind.closer <- which(d.cons.points < target.depth)
    Upperbounds[ind.closer, ] <- considered.points[ind.closer, ]
    considered.points[ind.closer, ] <- 0.5 * (Lowerbounds[ind.closer, ] +
                                            considered.points[ind.closer, ])
    #Find out which need to be further away from the center
    ind.further <- which(d.cons.points > target.depth)
    Lowerbounds[ind.further, ] <- considered.points[ind.further, ]
    considered.points[ind.further, ] <- 0.5 * (considered.points[ind.further,
                                                                  ] +
                                             Upperbounds[ind.further, ])

    #Add a breaking condition for extreme slow convergence situations
    temp1 <- sqrt(rowSums(matrix(abs(Lowerbounds - Upperbounds),
                                 ncol = NVar) ^ 2
                          )
                  )
    temp2 <- sqrt(rowSums(Upperbounds ^ 2))
    aid.ind <- which(temp1 / temp2 <= 10 ^ -10)
    ind.improve[aid.ind] <- FALSE

    count <- count + 1
  }

  considered.points <- sweep(x = considered.points,
                             MARGIN = 2, Center, FUN = "+")
  Converged[ind.improve] <- FALSE

  return(list(vertices = considered.points,
              converged = Converged)
  )
}
