
mrfDepth_theme <- function() {
  A <- theme_classic(base_size = 12, base_family = "") +
    theme(panel.background =  element_rect(fill = NA,
                                           colour = "black",
                                           size = 0.25),
          plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")
    )
  return(A)
}

GridPlot <- function(plotlist = NULL, layout = NULL) {
  num.plots <- length(plotlist)
  n.row <- nrow(layout)
  n.col <- ncol(layout)

  if (num.plots == 1) {
    print(plotlist[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(n.row, n.col)))

    # Make each plot, in the correct location
    for (i in 1:n.row) {
      for (j in 1:n.col) {
        if (layout[i, j] != 0)
        print(plotlist[[ layout[i, j] ]], vp = viewport(layout.pos.row = i,
                                                           layout.pos.col = j))
      }
    }
  }
}

matSubstract.c <- compiler:::cmpfun(function(mat, Center, n.row) {
  mat - rep(1, n.row) %*% Center
})


# pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 4, 4), "null"))))
# grid.text("title of this panel", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
# print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
# print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))