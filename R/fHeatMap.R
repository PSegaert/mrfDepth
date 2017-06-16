fheatmap <- function(rowValues, cellValues, type, scalename = ""){

  ######
  # Check input.
  if (missing(rowValues)) {
    stop("Input argument rowValues is required.")
  }
  if (is.null(rowValues)) {
    stop("Input argument rowValues is required.")
  }
  if (missing(cellValues)) {
    stop("Input argument cellValues is required.")
  }
  if (is.null(cellValues)) {
    stop("Input argument cellValues is required.")
  }
  if (missing(type)) {
    stop("Input argument type is required.")
  }


  x <- 1:ncol(cellValues)
  y <- nrow(cellValues):1
  xy <- expand.grid(x = x, y = y)

  if (type == "depth") {
    Ind <- sort(rowValues, decreasing = TRUE, index.return = TRUE)$ix
  }
  else{
    Ind <- sort(rowValues, decreasing = TRUE, index.return = TRUE)$ix
  }

  SortedcellValues <- cellValues[Ind, ]

  PlotData <- as.data.frame(cbind(xy, as.vector(t(SortedcellValues))))
  colnames(PlotData) <- c("x", "y", "Score")

  Score <- NULL
  Plot <- ggplot()
  Plot <- Plot + geom_tile(data = PlotData, aes(x, y, fill = Score))
  Plot <- Plot + coord_cartesian(xlim = c(-0.5, ncol(cellValues) + 1.5),
                                 ylim = c(-0.5, nrow(cellValues) + 1.5))
  if (type == "depth") {
    Plot <- Plot + scale_fill_gradientn(name = scalename,
                                       colours = c(low = "white",
                                                   high = "darkgreen"))
  }
  else{
    Plot <- Plot + scale_fill_gradientn(name = scalename,
                                       colours = c(low = "white",
                                                   high = "darkred"))
  }
  Plot <- Plot + scale_y_discrete(breaks = 1:length(Ind),
                                  labels = rev(Ind),
                                  name = "")
  Plot <- Plot + mrfDepth_theme() +
                 scale_x_continuous(name = "")

  return(Plot)
}
