#' Visualize matrix data in a form of a heatmap 
#' 
#' @description  
#' Produces convenient matrix data visualization in the form of a heatmap, with
#' the use of \code{ggplot2} library. 
#' 
#' @details 
#' Designed to provide 
#' decent visualization output with low number of input arguments that a user has to define.
#' Allows for further output plot adjustments that one often seeks for, including
#' adding a title, font size change, axis label clearing, axis label rotating etc.
#' 
#' By default, adjusts color scale based on matrix values: 
#' outputs 1-color plot for a matrix input which values are of one sign and
#' outputs 2-color plot for a matrix input which values are of both 
#' negative and positive values. 
#' 
#' By default, adjusts the plot legend color so as it is: 
#' white-red for a data matrix with only positive values,  
#' blue-white for a data matrix with only negative values,
#' blue-white-red for a data matrix with both positive and negative values.
#' 
#' By default, adjusts the color legend scale so as 
#' it starts from 0 for matrix values with one sign or is symmetric for matrix 
#' with both negative and positive values.
#' 
#' Does not return any object as it calls a \code{ggplot::plot} method itself.   
#' 
#' @param matrix.object data matrix 
#' @param title plot title 
#' @param base_size \code{ggplot2} theme \code{base_size} parameter value
#' @param limits 2-element vector of color legend scale limit values
#' 
#' @param colors.palette  color palette used as an argument in the 
#' \code{scale_fill_gradient} \code{ggplot2} expression
#' @param adjust.limits logical whether or not to adjust color legend scale limits 
#' in a way described in Details section
#' (if FALSE, default \code{ggplot2} options are applied)
#' @param adjust.colors logical whether or not to adjust plot color selection 
#' in a way described in Details section
#' (if FALSE, default \code{ggplot2} options are applied)
#' @param rotate.x.labels logical whether or not to rotate x-axis labels by 90 degrees
#' @param uniform.labes logical whether or not to add "generic" column and rows labeling
#'                        ("c1", "c2", ... for columns and "r1", "r2", ... for rows);
#'                        might be useful i.e. if the matrix we pass to function has 
#'                        some really long colnames and rownames assigned
#' @param clear.x.label logical whether or not to clear x-axis labels 
#' @param clear.y.label logical whether or not to clear y-axis labels 
#' @param clear.labels logical whether or not to clear both x- and y-axis labels     
#' @param geom_tile.colour \code{ggplot2} parameter defining the color of blocks grid 
#' 
#' @return \code{NULL}
#' 
#' @examples 
#' set.seed(123)
#' set.seed(123)
#' n <- 100
#' p <- 50
#' colnames <- replicate(p, paste0(sample(LETTERS, 5), collapse = ""))
#' rownames <- replicate(n, paste0(sample(LETTERS, 5), collapse = ""))
#' mat1 <- matrix(rnorm(n*p), nrow = n, ncol = p, dimnames = list(rownames, colnames))
#' vizu.mat(mat1, "mat1")
#' vizu.mat(mat1, "mat1 with uniformed labels", uniform.labes = TRUE, clear.labels = FALSE)
#' vizu.mat(mat1, "mat1 with predefined scale limits\nand smaller font size", 
#' limits = c(-4, 4), base_size = 5)
#' vizu.mat(cor(mat1), "column correlation matrix", geom_tile.colour = "white")
#' 
#' @import reshape2
#' @import ggplot2
#' @export
#' 
vizu.mat <- function(matrix.object, 
                     title = "", 
                     base_size = 12, 
                     limits = NULL, 
                     colors.palette = NULL,
                     adjust.limits = TRUE, 
                     adjust.colors = TRUE,
                     rotate.x.labels = TRUE, 
                     # change.pos.x.labels = TRUE,
                     uniform.labes = FALSE,
                     clear.x.label = FALSE, 
                     clear.y.label = FALSE, 
                     clear.labels = TRUE,
                     geom_tile.colour = "grey"){
                  
  # Correct colnames and rownames
  n.col <- ncol(matrix.object)
  n.row <- nrow(matrix.object)
  if (is.null(colnames(matrix.object)) | uniform.labes){
    colnames(matrix.object) <- paste0("c", 1:n.col)
  }
  if (is.null(rownames(matrix.object)) | uniform.labes){
    rownames(matrix.object) <- paste0("r", 1:n.row)
  }
  # Reshape matrix 
  matrix.object.m <- reshape2::melt(matrix.object)
  matrix.object.m[, "Var1"] <- factor(matrix.object.m[, "Var1"], levels = rev(rownames(matrix.object)))
  matrix.object.m[, "Var2"] <- factor(matrix.object.m[, "Var2"], levels = colnames(matrix.object))
  # Define plot base
  plot.tmp <- 
    ggplot(matrix.object.m, aes_string("Var2", "Var1")) + 
    geom_tile(ggplot2::aes_string(fill = "value"), colour = geom_tile.colour) + 
    labs(x = "", y = "", title = title) + 
    theme_grey(base_size = base_size)
  #     theme_grey(base_size = base_size) + 
  #     theme(panel.grid.minor = element_line(colour = "black", size = 5.5), 
  #           panel.grid.major = element_line(colour = "black", size = 5.2))
  # Adjust plot limits 
  mat.r <- range(matrix.object, na.rm = T)
  if (!is.null(limits)){
    limits.tmp <- limits
  } else if (adjust.limits){
    if (mat.r[1]<0 & mat.r[2]>0){
      limits.tmp <- c(-mat.r[2], mat.r[2])
      colours.tmp <- c("blue", "white", "red")
    } else if (mat.r[2] > 0){
      limits.tmp <- c(0, mat.r[2])
      colours.tmp <- c("white", "red")
    } else {
      limits.tmp <- c(mat.r[1], 0)
      colours.tmp <- c("blue", "white")
    }
  } else {
    limits.tmp <- NULL
  }
  # Adjust plot colors.palette 
  if (!is.null(colors.palette)){
    colours.tmp <- colors.palette
  } else if (adjust.colors){
    if (mat.r[1]<0 & mat.r[2]>0){
      colours.tmp <- c("blue", "white", "red")
    } else if (mat.r[2] > 0){
      colours.tmp <- c("white", "red")
    } else {
      colours.tmp <- c("blue", "white")
    }
  } else {
    colours.tmp <- NULL
  }
  plot.tmp <- plot.tmp + 
    scale_fill_gradientn(colours = colours.tmp, limits = limits.tmp)
  # Rotate x axis 
  if (rotate.x.labels){
    plot.tmp <- plot.tmp + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  # Clear axis labels and title 
  if (clear.labels){
    plot.tmp <- plot.tmp + 
      scale_x_discrete(breaks=NULL) + 
      scale_y_discrete(breaks=NULL)
  } else {
    if (clear.x.label) {
      plot.tmp <- plot.tmp + scale_x_discrete(breaks=NULL)
    }
    if (clear.y.label) {
      plot.tmp <- plot.tmp + scale_y_discrete(breaks=NULL)
    }
  } 
  # # Change x axis position to the top 
  # if (change.pos.x.labels & !clear.labels){
  #   plot.tmp <- ggdraw(switch_axis_position(plot.tmp, axis = 'x'))
  # }
  print(plot.tmp)
}