#' Size of Energy Trees
#' 
#' Depth and width of an Energy Tree.
#' 
#' @aliases depth width
#' 
#' @param x An Energy Tree object.
#' @param ... Additional arguments.
#' 
#' @name etree-size
NULL


#' @describeIn etree-size Depth of the three.
#' @export
depth <- function(x, ...)
  UseMethod("depth")


# @describeIn etree-size Number of terminal nodes in the tree.
#' @export
width <- function(x, ...)
  UseMethod("width")