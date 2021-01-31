

# A helper function for plotting results.
myplot <- function(object, title) {
  require(ggplot2)
  y_min <- 0.8*min(object$time)
  y_max <-  1.05 * max(object$time)
  object$ntime <- microbenchmark:::convert_to_unit(object$time, "t")
  plt <- ggplot(object, aes_string(x = "expr", y = "ntime", fill = "expr")) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    stat_ydensity() + scale_x_discrete(name = "")

  y_label <- sprintf("Time [%s]", attr(object$ntime,
                                       "unit"))

  plt <-  plt + scale_y_log10(name = y_label) + coord_flip() +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size =12) ) +
    labs(title = title) +
    theme(plot.title = element_text(size = 20, hjust = 0.5))
  plt
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
