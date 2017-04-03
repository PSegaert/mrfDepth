bagplot <- function(compbag.result,
                    colorbag = NULL,
                    colorloop = NULL,
                    colorchull = NULL,
                    databag = TRUE,
                    dataloop = TRUE,
                    plot.fence = FALSE){
  ######
  # Check input.
  if (missing(compbag.result)) {
    stop("Input argument compbag.result is required.")
  }

  if (!("mrfDepth" %in% class(compbag.result))) {
    stop("The input parameter compbag.result should be the
         return of a call to compBagplot. ")
  }
  if (!("compBagplot" %in% class(compbag.result))) {
    stop("The input parameter compbag.result should be the
         return of a call to compBagplot. ")
  }
  if (is.null(colorbag)) {
    colorbag <- rgb(0.6, 0.6, 1)
  }
  if (is.null(colorloop)) {
    colorloop <- rgb(0.8, 0.8, 1)
  }
  if (is.null(colorchull)) {
    colorchull <- rgb(1, 1, 1)
  }

  ind.bag <- which(compbag.result$datatype[, 3] == 1)
  data.bagcontour <- data.frame(compbag.result$bag)
  data.inbag <- data.frame(compbag.result$datatype[ind.bag, 1:2,
                                                    drop = FALSE])
  ind.loop <- which(compbag.result$datatype[, 3] == 2)
  data.inloop <- data.frame(compbag.result$datatype[ind.loop, 1:2,
                                                     drop = FALSE])
  data.infence <- data.frame(compbag.result$fence)
  ind.outl <- which(compbag.result$datatype[, 3] == 3)
  data.outliers <- data.frame(compbag.result$datatype[ind.outl, 1:2,
                                                      drop = FALSE])
  colnames(data.infence) <- c("x", "y")
  colnames(data.inloop) <- c("x", "y")
  colnames(data.inbag) <- c("x", "y")
  colnames(data.bagcontour) <- c("x", "y")
  colnames(data.outliers) <- c("x", "y")

  label.x <- colnames(compbag.result$datatype)[1]
  label.y <- colnames(compbag.result$datatype)[2]

  # Initialise plot
  plot <- ggplot()

  # Add the fence
  if (plot.fence) {
    data.reg <- data.infence[chull(data.infence), ]
    colnames(data.reg) <- c("x", "y")
    data.reg <- rbind(data.reg, data.reg[1, ])
    plot <- plot + geom_path(data = data.reg,
                             mapping = aes_string(x = "x", y = "y"),
                             linetype = "dashed"
                             )
  }

  # Add the filled loop
  data.reg <- data.inloop[chull(data.inloop), ]
  colnames(data.reg) <- c("x", "y")
  plot <- plot + geom_polygon(data = data.reg,
                             mapping = aes_string(x = "x", y = "y"),
                             fill = colorloop)

  # Add the filled bag
  data.reg <- data.bagcontour[chull(data.bagcontour), ]
  colnames(data.reg) <- c("x", "y")
  plot <- plot + geom_polygon(data = data.reg,
                             mapping = aes_string(x = "x", y = "y"),
                             fill = colorbag)

  # Add the innermost contour if applicable
  if (nrow(compbag.result$chull) > 1) {
    data.chull <- data.frame(compbag.result$chull)
    colnames(data.chull) <- c("x", "y")
    plot <- plot + geom_polygon(data = data.chull,
                               mapping = aes_string(x = "x", y = "y"),
                               fill = colorchull)
  }

  # Add the center
  data.reg <- data.frame(matrix(compbag.result$center, ncol = 2))
  colnames(data.reg) <- c("x", "y")
  plot <- plot + geom_point(data = data.reg,
                           mapping = aes_string(x = "x", y = "y"),
                           shape = 23, size = 4, color = "red", fill = "red"
                           )
  # If requested plot the points
  if (databag) {
    plot <- plot + geom_point(data = data.inbag,
                             mapping = aes_string(x = "x", y = "y")
                    )
  }
  if (dataloop) {
    plot <- plot + geom_point(data = data.inloop,
                             mapping = aes_string(x = "x", y = "y")
                             )
  }

  # Plot outliers
  if (nrow(data.outliers) > 0) {
    plot <- plot + geom_point(data = data.outliers,
                             mapping = aes_string(x = "x", y = "y"),
                             shape = 8, color = "red"
                             )
  }

  # Set up the figure
  plot.data <- compbag.result$datatype[, -3]
  colnames(plot.data) <- c("x", "y")
  plot.data <- rbind(plot.data, data.infence)
  x.range <- extendrange(plot.data[, 1], f = 0.05)
  y.range <- extendrange(plot.data[, 2], f = 0.05)
  plot <- plot + coord_cartesian(xlim = x.range, ylim = y.range)

  # give plot the package look
  plot <- plot + mrfDepth_theme()

  # Finalise
  plot <- plot + xlab(label.x) + ylab(label.y)
  if (compbag.result$type == "hdepth") {
    label.title <- paste("Bagplot based on halfspace depth")
  } else if (compbag.result$type == "projdepth") {
    label.title <- paste("Bagplot based on projection depth")
  } else if (compbag.result$type == "sprojdepth") {
    label.title <- paste("Bagplot based on skewness-adjusted projection depth")
  } else{
    label.title <- paste("Bagplot based on", compbag.result$type, "depth")
  }
  plot <- plot + ggtitle(label.title)

  return(plot)

}
