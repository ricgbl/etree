#' Size of Energy Trees
#' 
#' Depth and width of an Energy Tree.
#' 
#' @aliases depth width
#' 
#' @param x An object of class \code{etree}.
#' @param ... Additional arguments.
#' 
#' @returns
#' \code{depth()} returns the depth of the tree and \code{width()} gives the 
#' number of terminal nodes.
#' 
#' @name etree-size
NULL


#' @describeIn etree-size Depth of the three.
#' @export
depth <- function(x, ...)
  UseMethod("depth")


#' @describeIn etree-size Number of terminal nodes in the tree.
#' @export
width <- function(x, ...)
  UseMethod("width")