plotContours <- function(x, depthContour, data = TRUE) {

  ######
  # Check input.
  if (missing(x)) {
    stop("The input argument x is required.")
  }
  if (missing(depthContour)) {
    stop("The input argument depthContour is required.")
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
    stop("The data matrix x must be two dimensional.")
  }
  if (is.null(colnames(x))) {
    colnames(x) <- c("variable 1", "variable 2")
  }

  if (!("mrfDepth" %in% class(depthContour))) {
    stop("depthContour should be the return of a call to depthContour.")
  }
  if (!("depthContour" %in% class(depthContour))) {
    stop("depthContour should be the return of a call to depthContour.")
  }

  # Initialise plot
  plot <- ggplot()

  # If requested plot the points
  x.data <- data.frame(x)
  colnames(x.data) <- c("x", "y")
  if (data) {
    plot <- plot + geom_point(data = x.data,
                              mapping = aes_string(x = "x", y = "y")
                              )
  }

  # Add the contours
  for (i in 1:(length(depthContour) - 1)) {
    TResult <- depthContour[[i]]
    if (TResult$empty != 1) {
      # Add the filled bag
      data.cont <- data.frame(TResult$vertices[chull(TResult$vertices), ])
      data.cont <- rbind(data.cont, data.cont[1, ])
      colnames(data.cont) <- c("x", "y")
      plot <- plot + geom_path(data = data.cont,
                                  mapping = aes_string(x = "x", y = "y")
                                  )
    }
  }

  # give plot the package look
  plot <- plot + mrfDepth_theme()

  # Finalise
  plot <- plot + xlab(colnames(x)[1]) + ylab(colnames(x)[2])
  if (depthContour$type == "hdepth") {
    title.label <- paste("Halfspace depth contours")
  } else if (depthContour$type == "projdepth") {
    title.label <- paste("Projection depth contours")
  } else if (depthContour$type == "sprojdepth") {
    title.label <- paste("Skew-adjusted projection depth contours")
  } else{
    title.label <- ""
  }
  plot <- plot + ggtitle(title.label)


  return(plot)

}
