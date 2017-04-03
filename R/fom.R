fom <- function(fOutlResult, cutoff = FALSE){
  ######
  # Check input.
  if (missing(fOutlResult)) {
    stop("Input argument fOutlResult is required.")
  }

  #Check fOutlResult
  if (!is.list(fOutlResult)) {
    stop("fOutlResult must be a list returned from a call to fOutl.")
  }
  InputNames <- names(fOutlResult)
  if (!("fOutl" %in% class(fOutlResult))) {
    stop("fOutlResult must be a list returned from a call to fOutl.")
  }
  if (!("distType" %in% InputNames)) {
    stop("fOutlResult must be a list returned from a call to fOutl.")
  }
  if (!("weights" %in% InputNames)) {
    stop("fOutlResult must be a list returned from a call to fOutl.")
  }
  if (!("crossDistsX" %in% InputNames)) {
    stop(paste("fOutlResult must be a list returned from a call to fOutl",
               "with option diagnostic = TRUE"))
  }
  if (!("locOutlX" %in% InputNames)) {
    stop(paste("fOutlResult must be a list returned from a call to fOutl",
               "with option diagnostic = TRUE"))
  }

  #Get diagnostic info
  NFunc <- nrow(fOutlResult$crossDistsX)
  NTObs <- ncol(fOutlResult$crossDistsX)

  #Make local variables
  AOValues <- fOutlResult$crossDistsX
  fAO <- AOValues %*% fOutlResult$weights
  LocOutl <- rowSums(fOutlResult$locOutlX)

  PlotData <- data.frame(row.names = 1:NFunc)
  PlotData$fAO <- fAO
  PlotData$DispMeasure <- rowSds(AOValues) / (1 + fAO)
  PlotData$Shape <- rep(16, NFunc)
  PlotData$Shape[which(LocOutl >= 0.25 * NTObs)] <- 17
  PlotData$Shape[which(LocOutl >= 0.5 * NTObs)] <- 15
  PlotData$Shape[which(LocOutl >= 0.75 * NTObs)] <- 18
  PlotData$colorvec <- rep("black", NFunc)

  if (cutoff) {
    # Construct Cutoff
    CAO <- log(0.1 + sqrt(
                      (PlotData$fAO / median(PlotData$fAO)) ^ 2 +
                      (PlotData$DispMeasure / median(PlotData$DispMeasure)) ^ 2)
              )
    Fence <- qnorm(0.995) * mad(CAO) + median(CAO)
    theta <- seq(0, pi / 2, length = (100))
    FenceData <- matrix(0.0, nrow = length(theta), ncol = 2)
    colnames(FenceData) <- c("x", "y")
    FenceData <- data.frame(FenceData)
    FenceData$x <- median(PlotData$fAO) * (exp(Fence) - 0.1) * cos(theta)
    FenceData$y <- median(PlotData$DispMeasure) * (exp(Fence) - 0.1) * sin(theta)

    # Add colors
    PlotData$colorvec <- rep("black", length(CAO))
    PlotData$colorvec[which(CAO > Fence)] <- "red"
  }

  # Produce actual plot
  Plot <- ggplot(PlotData) + scale_shape_identity()
  Plot <- Plot + geom_point(mapping = aes_string(x = "fAO", y = "DispMeasure"),
                            color = PlotData$colorvec,
                            shape = PlotData$Shape,
                            size = 2) 
  if (cutoff) {
    Plot <- Plot + geom_path(mapping = aes_string(x = "x", y = "y"),
                             data = FenceData,
                             color = "black",
                             linetype = 2,
                             size = 1)
  }
  Plot <- Plot + xlab(fOutlResult$distType)
  Plot <- Plot + ylab(paste("sd(", fOutlResult$distType,
                            ") / (1+", fOutlResult$distType, ")",
                            sep = ""))
  Plot <- Plot + mrfDepth_theme() + guides(shape = guide_legend(title = "z"))
  Plot

  return(Plot)
}
