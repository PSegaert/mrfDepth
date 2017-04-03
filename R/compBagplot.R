compBagplot <- function(x, type="hdepth", sizesubset=500,
                        extra.directions = FALSE, options=NULL){

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
  n <- nrow(x)
  p <- ncol(x)
  if (n > sum(complete.cases(x))) {
    stop("Missing values in x are not allowed.")
  }
  if (p != 2) {
    stop("The bagplot can only be drawn for bivariate data.")
  }
  if (n < 10) {
    stop("At least 10 data points are required.")
  }
  if (is.null(colnames(x))) {
    colnames(x) <- c("variable 1", "variable 2")
  }

  # Check type
  Indtype <- match(type, c("hdepth", "projdepth", "sprojdepth"))[1]
  if (is.na(Indtype)) {
    stop("The input parameter type should be one of hdepth,
         projdepth or sprojdepth.")
  }

  # Check sizesubset
  if (!is.numeric(sizesubset)) {
    stop("The input parameter sizesubset should be numeric.")
  }

  # Check extra.directions
  if (!is.logical(extra.directions)) {
    stop("The input parameter extra.directions should be a logical.")
  }


  # check options
  if (is.null(options)) {
    options <- list()
  }
  if (!is.list(options)) {
    stop("options must be a list")
  }
  if ("max.iter"  %in% names(options)) {
    max.iter <- options[["max.iter"]]
    if (is.numeric(max.iter)) {
      if (max.iter < 1) {
        stop("The option parameter max.iter must be a strictly
             positive integer.")
      }
    }
    if (!is.numeric(max.iter)) {
      stop("The option parameter max.iter must be a strictly
           positive integer.")
    }
  } else {
    options$max.iter <- 100
  }


  #####
  # Check data for possible exact fit situations.
  tol <- 1e-7
  scaled.x <- scale(x)
  temp <- attributes(scaled.x)
  column.sd <- temp[["scaled:scale"]]
  if (sum(column.sd <= 1e-14) > 0) {
    warning("At least one of the variables has zero standard deviation.
            Check the data matrix x.")
    returned.result <- list(center = NULL,
                            chull = NULL,
                            bag = NULL,
                            fence = NULL,
                            datatype = NULL,
                            flag = NULL,
                            depth = NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14),
                            type = type)
    class(returned.result) <- c("mrfDepth", "compBagplot")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x
            could not be computed.")
    returned.result <- list(center = NULL,
                            chull = NULL,
                            bag = NULL,
                            fence = NULL,
                            datatype = NULL,
                            flag = NULL,
                            depth = NULL,
                            dimension = NULL,
                            hyperplane = NULL,
                            type = type)
    class(returned.result) <- c("mrfDepth", "compBagplot")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more details.")
    returned.result <- list(center = NULL,
                            chull = NULL,
                            bag = NULL,
                            fence = NULL,
                            datatype = NULL,
                            flag = NULL,
                            depth = NULL,
                            dimension = sum(w1$d > tol),
                            hyperplane = w1$v[, which(w1$d == min(w1$d))[1]],
                            type = type)
    class(returned.result) <- c("mrfDepth", "compBagplot")
    return(returned.result)
  }

  if (type == "hdepth") {
    result <- CompHalfSpaceBagplot(x, sizesubset)
  }
  else if (type == "projdepth" || type == "sprojdepth") {
    result <- CompProjectionTypeBagplot(x, type, extra.directions, options)
  }
  else{
    stop("The input parameter type is not recognized.")
  }

  class(result) <- c("mrfDepth", "compBagplot")
  
  return(result)

}

CompHalfSpaceBagplot <- function(x, sizesubset) {
  n <- nrow(x)
  interpol <- matrix(0.0, 2 * n, 2)
  datatyp <- matrix(0.0, n, 3)
  datatyp2 <- matrix(0.0, n, 2)
  chull <- matrix(0.0, sizesubset * (sizesubset - 1) / 2, 2)
  pxpy <- matrix(0.0, n, 3)
  PDepths <- matrix(0.0, n, 2)
  storage.mode(interpol) <- "double"
  storage.mode(datatyp) <- "double"
  storage.mode(datatyp2) <- "double"
  storage.mode(chull) <- "double"
  storage.mode(pxpy) <- "double"
  storage.mode(PDepths) <- "double"

  result <- .Fortran("bagplotf",
                     as.integer(n),          #1 Number of points
                     as.double(x[, 1]),      #2 First variable of x
                     as.double(x[, 2]),      #3 Second variable of x
                     as.integer(3),          #4 Specify whisker format
                     as.integer(sizesubset), #5 nsub
                     as.double(rep(0, 2)),   #6 tukm
                     as.integer(1),          #7 Number of points in chull
                     chull,                  #8 Contour points of deepest region
                     interpol,               #9 Points in the bag
                     as.integer(0),          #10  Number of points in bag
                     datatyp,                #11  Data with their type
                     as.integer(rep(0, n)),  #12  Indices of outliers
                     as.integer(0),          #13  Number of outliers
                     datatyp2,               #14  datatyp2
                     pxpy,                   #15  Aid variable
                     as.integer(1),          #16  Flag indicating if half of the
                                             #    points lie on a vertical line
                     as.integer(1),          #17 Flag to indicate whether
                                             #   dithering was performed
                     PDepths,                #18 Depths of points
                     PACKAGE = "mrfDepth")

  if (result[[13]] > 0) {
     ind.outl <- result[[12]][1:result[[13]]]
  } else {
    ind.outl <- NULL
  }
  Flag <- rep(0, n)
  Flag[ind.outl] <- 1

  if (result[[16]] == 1) {
       stop("At least half of the observations lie on a line,
            a univariate boxplot should be used.")
  }

  Datatyp <- result[[11]]
  colnames(Datatyp) <- c(colnames(x), "positionIndicator")

  # Calculate the fence.
  Bag <- result[[9]][1:result[[10]], ]
  Center <- result[[6]]
  centered.bag <- sweep(Bag, MARGIN = 2, Center, FUN = "-")
  fence <-  3 * centered.bag
  fence <- sweep(fence, MARGIN = 2, Center, FUN = "+")

  returned.result <- list(center = result[[6]],
                chull = result[[8]][1:result[[7]], ],
                bag = result[[9]][1:result[[10]], ],
                fence = fence,
                datatype = Datatyp,
                flag = Flag,
                depth = result[[18]][, 1],
                dimension = NULL,
                hyperplane = NULL,
                type = "hdepth"
                )
  class(returned.result) <- c("mrfDepth", "compBagplot")
  return(returned.result)
}

CompProjectionTypeBagplot <- function(x, type, extra.directions, options){

  if (type == "projdepth") {
    Result <- projdepth(x = x, z = x, options = options)
    Center <- projmedian(x = x, projection.depths = Result[["depthZ"]])$max
  }
  else if (type == "sprojdepth") {
    Result <- sprojdepth(x = x, z = x, options = options)
    Center <- sprojmedian(x = x, sprojection.depths = Result[["depthZ"]])$max
  }
  else{
    stop("The input parameter type is not recognized.")
  }

  n <- nrow(x)
  p <- ncol(x)
  n.bag <- ceiling(n / 2)
  ind.outl <- (Result[["depthZ"]] < Result$cutoff)

  sorted.depth <- sort(Result[["depthZ"]], decreasing = TRUE)
  bag.cutoff <- sorted.depth[n.bag]
  ind.bag <- which(Result[["depthZ"]] >= bag.cutoff)
  Data1 <- x[ind.bag, 1:2]
  TInd <- chull(Data1[, 1], Data1[, 2])
  Bag <- Data1[TInd, ]

  Datatyp <- cbind(x, rep(0.0, n))
  Datatyp[!ind.outl, 3] <- 2
  Datatyp[ind.bag, 3] <-  1
  Datatyp[ind.outl, 3] <- 3

  colnames(Datatyp) <- c(colnames(x), "positionIndicator")

  TInd <- which(Datatyp[, 3] != 3)
  ind.loop <- chull(x = Datatyp[TInd, 1], y = Datatyp[TInd, 2])
  Loop <- Datatyp[TInd[ind.loop], 1:2]
  Directions <- sweep(rbind(Loop, Bag), MARGIN = 2, Center, FUN = "-")
  if (extra.directions) {
    Angle <- seq(0, 2 * pi, length.out = 250)
    Directions <- rbind(Directions, cbind(cos(Angle), sin(Angle)))
  }
  Directions <- Directions / sqrt(rowSums(Directions ^ 2))

  fence <- CalculateOneProjTypeContour(x = x,
                                      Center = Center,
                                      alpha = Result$cutoff,
                                      type = type,
                                      Directions = Directions,
                                      options = options)
  fence <- fence$vertices

  result <- list(center = matrix(Center, ncol = p),
                 chull = matrix(Center, ncol = p),
                 bag = Bag,
                 fence = fence,
                 datatype = Datatyp,
                 flag = as.numeric(Result[["flagZ"]]),
                 depth = Result[["depthZ"]],
                 dimension = NULL,
                 hyperplane = NULL,
                 type = type)

  class(result) <- c("mrfDepth", "compBagplot")
  return(result)

}
