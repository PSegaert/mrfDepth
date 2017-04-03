bagdistance <- function(x, z = NULL, options = list()){

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
        stop("Option max.iter must be a strictly positive integer.")
      }
    }
    if (!is.numeric(max.iter)) {
      stop("Option max.iter must be a strictly positive integer.")
    }
  } else {
    options$max.iter <- 100
  }
  if ("approx"  %in% names(options)) {
    approx <- options[["approx"]]
    if (!is.logical(approx)) {
      stop("Option approx must be a logical value.")
    }
  } else {
    options$approx <- TRUE
  }

  #####
  # Check data for possible exact fit situations.
  tol <- 1e-7
  scaled.x <- scale(x)
  temp <- attributes(scaled.x)
  column.sd <- temp[["scaled:scale"]]
  if (sum(column.sd <= 1e-14) > 0) {
    warning("One of the variables of x has zero standard deviation.
            Check the data matrix x.")
    returned.result <- list(bagdistance = NULL,
                            cutoff = NULL,
                            flag = NULL,
                            converged = NULL,
                            dimension = sum(column.sd > 1e-14),
                            hyperplane = as.numeric(column.sd <= 1e-14))
    class(returned.result) <- c("mrfDepth", "bagdistance")
    return(returned.result)
  }
  w1 <- try(svd(scaled.x / sqrt(n1 - 1)), silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x
            could not be computed.")
    returned.result <- list(bagdistance = NULL,
                            cutoff = NULL,
                            flag = NULL,
                            converged = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "bagdistance")
    return(returned.result)
  }
  if (min(w1$d) < tol) {
    warning("An exact fit is found. Check the output for more details.")
    returned.result <- list(bagdistance = NULL,
                            cutoff = NULL,
                            flag = NULL,
                            converged = NULL,
                            dimension = sum(w1$d > tol),
                            hyperplane = w1$v[, which(w1$d == min(w1$d))[1]])
    class(returned.result) <- c("mrfDepth", "bagdistance")
    return(returned.result)
  }
  #####
  # Prepare the start of the algorithm
  #####
  # Find the median of the depths of the data points.
  result1 <- hdepth(x = x, z = x, options = options)
  if (is.null(result1[["depthZ"]])) {
    stop("The halfspace depth of x can not be computed.")
  }
  # Find the halfspace median of x.
  result2 <- hdepthmedian(x = x)
  if (is.null(result2[["median"]])) {
    stop("The halfspace median of x can not be computed.")
  }
  center <- result2[["median"]]

  # Initialise constants.
  n.var <- p1
  # Note that in the hdepth call z was taken to be x.
  target.depth <- sort(result1[["depthZ"]], decreasing = TRUE)[ceiling(n1 / 2)]
  n.points <- n2

  # Initialise vector of halfspace distances.
  distances <- rep(0.0, n.points)

  #####
  # Start the real calculations.
  #####
  if (n.var == 1) {
    # This case is pretty easy.
    depths <- result1[["depthX"]]
    Ind <- which(depths >= target.depth)
    lower.bound <- min(x[Ind])
    upper.bound <- max(x[Ind])

    if (abs(lower.bound - center) < tol | abs(upper.bound - center) < tol) {
      warning("The bag has length zero.")
      returned.result <- list(bagdistance = NULL,
                              cutoff = NULL,
                              flag = NULL,
                              converged = NULL,
                              dimension = NULL,
                              hyperplane = NULL)
      class(returned.result) <- c("mrfDepth", "bagdistance")
      return(returned.result)
    }

    Ind <- which(z <= center)
    distances[Ind] <- abs(z[Ind] - center) / abs(lower.bound - center)

    Ind <- which(z >= center)
    distances[Ind] <- abs(z[Ind] - center) / abs(upper.bound - center)
    cutoff <- sqrt(qchisq(0.99, p1))
    flag.z <- distances <= cutoff

    returned.result <- list(bagdistance = distances,
                cutoff = cutoff,
                flag = flag.z,
                converged = NULL,
                dimension = NULL,
                hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "bagdistance")
    return(returned.result)
  }
  else if (n.var == 2 & options$approx == FALSE) {
    # First calculate the isodepth contour and then calculate
    # the intersection of each ray with the contour.

    # Center all data
    data.distribution <- sweep(x, MARGIN = 2, center, FUN = "-")
    data.to.calc.dist <- sweep(z, MARGIN = 2, center, FUN = "-")

    # Compute the isodepth contour and calculate intersections.
    vertices.countour <- depthContour(data.distribution,
                                      alpha = target.depth,
                                      type = "hdepth"
                                      )$Contour
    if ("dataDimension" %in% names(vertices.countour)) {
      # An exact fit was found
      stop("The isodepth contour could not be computed.")
    }
    vertices.countour <- vertices.countour$vertices
    n.vertices <- nrow(vertices.countour)

    # Define an angle transformation function such that all points
    # have increasing angles in [0, 2pi)
      AngleTransfo <- function(angles, data){
      ind.180 <- which(data[, 1] < 0)
      ind.360 <- which(data[, 2] < 0 & data[, 1] > 0)

      angles[ind.180] <- angles[ind.180] + pi
      angles[ind.360] <- angles[ind.360] + (2 * pi)
      return(angles)
    }

    # Find the angels of all points.
    angles.contour <- atan(vertices.countour[, 2] / vertices.countour[, 1])
    angles.contour <- AngleTransfo(angles = angles.contour,
                                   data = vertices.countour)

    angles.distri.data <- atan(data.distribution[, 2] / data.distribution[, 1])
    angles.distri.data <- AngleTransfo(angles = angles.distri.data,
                                       data = data.distribution)

    angles.to.calc.dist <- atan(data.to.calc.dist[, 2] / data.to.calc.dist[, 1])
    angles.to.calc.dist <- AngleTransfo(angles = angles.to.calc.dist,
                                        data = data.to.calc.dist)

    # Sort all angles keeping track of which group the angles belong to.
    temp.angles <- rbind(cbind(angles.contour,
                               rep(1.0, length(angles.contour))),
                         cbind(angles.to.calc.dist,
                               rep(0.0, n.points))
                         )
    temp.data <- rbind(vertices.countour, data.to.calc.dist)
    temp.ind <- order(temp.angles[, 1], temp.angles[, 2], decreasing = FALSE)
    temp.angles <- temp.angles[temp.ind, ]
    temp.data <- temp.data[temp.ind, ]

    # Start calculating the intersections.
    ind.distri <- which(temp.angles[, 2] == 1)
    n.ind.distri <- length(ind.distri)

    # Do calculations for points between 2pi and 0.
    to.calc.ind <- c()
    if (ind.distri[1] != 1) {
      to.calc.ind <- 1:(ind.distri[1] - 1)
    }
    if (ind.distri[n.ind.distri] != length(temp.ind)) {
      to.calc.ind <- c(to.calc.ind,
                       (ind.distri[n.ind.distri] + 1):length(temp.ind)
                       )
    }
    if (!is.null(to.calc.ind)) {
      # Find the equation of the line element in the contour.
      if (temp.data[ind.distri[1], 1] ==
          temp.data[ind.distri[n.ind.distri], 1]) {
        # The line segment is vertical
        Temp <- temp.data[to.calc.ind, ]
        if (is.vector(Temp)) {
          Temp <- matrix(Temp, ncol = 2)
        }
        Intersections <- cbind(rep(temp.data[ind.distri[1], 1],
                                   ind.distri[1] - 1),
                              (Temp[, 2] / Temp[, 1]) *
                                temp.data[ind.distri[1], 1]
                              )
        distances[temp.ind[to.calc.ind] - n.vertices] <-
          sqrt(Temp[, 1] ^ 2 + Temp[, 2] ^ 2) /
          (sqrt(Intersections[, 1] ^ 2 + Intersections[, 2] ^ 2))
      }
      else{
        Temp <- temp.data[to.calc.ind, ]
        if (is.vector(Temp)) {
          Temp <- matrix(Temp, ncol = 2)
        }
        segment.slope <- (temp.data[ind.distri[1], 2] -
                           temp.data[ind.distri[n.ind.distri], 2]) /
                        (temp.data[ind.distri[1], 1] -
                           temp.data[ind.distri[n.ind.distri], 1])
        intersect.x <- (segment.slope * temp.data[ind.distri[1], 1] -
                         temp.data[ind.distri[1], 2]) /
                      (segment.slope - Temp[, 2] / Temp[, 1])
        intersect.y <- (Temp[, 2] / Temp[, 1]) * intersect.x
        Intersections <- cbind(intersect.x, intersect.y)
        distances[temp.ind[to.calc.ind] - n.vertices] <-
                          sqrt(Temp[, 1] ^ 2 + Temp[, 2] ^ 2) /
                          sqrt(Intersections[, 1] ^ 2 + Intersections[, 2] ^ 2)
      }
    }
    for (i in 1:(length(ind.distri) - 1)) {
      if ((ind.distri[i] + 1) != ind.distri[i + 1]) {
        to.calc.ind <- (ind.distri[i] + 1):(ind.distri[i + 1] - 1)
        # Find the equation of the line element in the contour.
        if (temp.data[ind.distri[i], 1] == temp.data[ind.distri[i + 1], 1]) {
          # The line segment is vertical
          Temp <- temp.data[to.calc.ind, ]
          if (is.vector(Temp)) {
            Temp <- matrix(Temp, ncol = 2)
          }
          Intersections <- cbind(rep(temp.data[ind.distri[i], 1], nrow(Temp)),
                                (Temp[, 2] / Temp[, 1]) *
                                  temp.data[ind.distri[i], 1]
                                )
          distances[temp.ind[to.calc.ind] - n.vertices] <-
            sqrt(Temp[, 1] ^ 2 + Temp[, 2] ^ 2) /
            sqrt(Intersections[, 1] ^ 2 + Intersections[, 2] ^ 2)
        }
        else{
          Temp <- temp.data[to.calc.ind, ]
          if (is.vector(Temp)) {
            Temp <- matrix(Temp, ncol = 2)
          }
          segment.slope <- (temp.data[ind.distri[i], 2] -
                              temp.data[ind.distri[i + 1], 2]) /
                           (temp.data[ind.distri[i], 1] -
                              temp.data[ind.distri[i + 1], 1])
          intersect.x <- (segment.slope * temp.data[ind.distri[i], 1] -
                            temp.data[ind.distri[i], 2]) /
                         (segment.slope - Temp[, 2] / Temp[, 1])
          intersect.y <- (Temp[, 2] / Temp[, 1]) * intersect.x
          Intersections <- cbind(intersect.x, intersect.y)
          dist.temp <- sqrt(Temp[, 1] ^ 2 + Temp[, 2] ^ 2) /
                      sqrt(Intersections[, 1] ^ 2 + Intersections[, 2] ^ 2)
          distances[temp.ind[to.calc.ind] - n.vertices] <- dist.temp
        }
      }
    }
    distances[is.nan(distances)] <- 0

    cutoff <- sqrt(qchisq(0.99, p1))
    flag.z <- distances <= cutoff

    returned.result <- list(bagdistance = distances,
                            cutoff = cutoff,
                            flag = flag.z,
                            converged = NULL,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "bagdistance")
    return(returned.result)
  }
  else{
    # Center all data
    data.distribution <- sweep(x, MARGIN = 2, center, FUN = "-")
    data.to.calc.dist <- sweep(z, MARGIN = 2, center, FUN = "-")

    # Get the directions
    directions <- data.to.calc.dist
    directions <-  directions / sqrt(rowSums(directions ^ 2))

    result <-  CalcOneHalfContourIntersect(x = data.distribution,
                                           Center = rep(0.0, n.var),
                                           alpha = target.depth,
                                           Directions = directions,
                                           options = options
                                           )
    intersects <- result$vertices

    distances <- sqrt(rowSums(data.to.calc.dist ^ 2))
    scales <- sqrt(rowSums(intersects ^ 2))
    scaled.distances <- distances / scales
    scaled.distances[is.nan(scaled.distances)] <- 0

    cutoff <- sqrt(qchisq(0.99, p1))
    flag.z <- scaled.distances <= cutoff

    returned.result <- list(bagdistance = scaled.distances,
                            cutoff = cutoff,
                            flag = flag.z,
                            converged = result$converged,
                            dimension = NULL,
                            hyperplane = NULL)
    class(returned.result) <- c("mrfDepth", "bagdistance")
    return(returned.result)
    }
}
