mrainbowplot <- function(x, depths, col = NULL, plot.options = list()) {

  ######
  # Check input.
  if (missing(x)) {
    stop("The input argument x is required.")
  }
  if (missing(depths)) {
    stop("The input argument depths is required.")
  }

  # Check x
  x <- data.matrix(x)
  if (!is.matrix(x)) {
    stop("The input argument x must be a matrix.")
  }
  if (sum(is.nan(x)) != 0) {
    stop("The data matrix x contains missing cases.")
  }
  n1 <- dim(x)[1]
  p1 <- dim(x)[2]
  if (p1 != 2) {
    stop("The data matrix x needs to be two-dimensional.")
  }
  plot.data <- data.frame(cbind(x, depths))
  colnames(plot.data) <- c("x", "y", "depth")
  if (is.null(colnames(x))) {
    Labs <- c("variable 1", "variable 2")
  } else {
    Labs <- colnames(x)
  }

  # Check depths
  depths <- data.matrix(depths)
  if (sum(is.nan(depths)) != 0) {
    stop("The input argument depths contains missing cases.")
  }
  n2 <- dim(depths)[1]
  p2 <- dim(depths)[2]
  if (n1 != n2) {
    stop("The input arguments depths and x do not have the same size.")
  }
  if (p2 != 1) {
    stop("The input argument depths must be column matrix.")
  }

  # Check rgb colors
  if (is.null(col)) {
    col <- makeColors_MRainbow()
  }
  col <- data.matrix(col)
  if (sum(is.nan(col)) != 0) {
    stop("The input argument col contains missing cases.")
  }
  p3 <- dim(col)[2]
  if (p3 != 3) {
    stop("The input argument col must have three columns.")
  }
  if (sum(col > 1) != 0) {
    stop("All values in the paramter col must lie in [0,1].")
  }
  if (sum(col < 0) != 0) {
    stop("All values in the paramter col must lie in [0,1].")
  }
  RGBCols <- rgb(col[, 1], col[, 2], col[, 3])

  #check plot.options and load defaults if necessary.
  if (!is.list(plot.options)) {
    stop("options must be a list")
  }
  if ("legend.title" %in% names(options)) {
    legend.title <- options[["legend.title"]]
  } else {
    legend.title <- "Depth"
  }
  if ("point.size" %in% names(options)) {
    point.size <- options[["point.size"]]
  } else {
    point.size <- 4
  }
  if(!is.character(legend.title)){
    stop("legend.title must be a string of characters.")
  }
  if(!is.numeric(point.size)){
    stop("point.size must be numeric.")
  } else{
    if(point.size <= 0){
      stop("point.size must be a strictly positive numeric.")
    }
  }
  

  # Create basic plot
  plot <- ggplot()
  plot <- plot + geom_point(data = plot.data,
                          mapping = aes_string(x = "x", y = "y",
                                                  colour = "depth"),
                          size = point.size
                             )
  # Make a colorbal scale
  plot <- plot + scale_colour_gradientn(colours = RGBCols, name = legend.title)
  plot

  # Set up the figure
  plot <- plot + ggtitle("Bivariate rainbow plot")
  x.range <- extendrange(x[, 1], f = 0.05)
  y.range <- extendrange(x[, 2], f = 0.05)
  plot <- plot + coord_cartesian(xlim = x.range, ylim = y.range)

  # give plot the package look
  plot <- plot + mrfDepth_theme() +
    scale_x_continuous(name = Labs[1]) +
    scale_y_continuous(name = Labs[2])

  return(plot)

}

makeColors_MRainbow <- function() {
  RGBmatrix <- c(0.823529411764706, 1.000000000000000, 0.823529411764706,
                 0.392156862745098, 0.823529411764706, 0.392156862745098,
                 0.235294117647059, 0.784313725490196, 0.235294117647059,
                 0.117647058823529, 0.666666666666667, 0.117647058823529,
                 0.058823529411765, 0.509803921568627, 0.058823529411765,
                                 0, 0.196078431372549, 0)
  RGBmatrix <- matrix(RGBmatrix, ncol = 3, byrow = TRUE)
}
