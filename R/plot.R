#' Visualization of Energy Trees
#'
#' Returns the plot of an object of class \code{"etree"}.
#'
#' @param x	An object of class \code{"etree"}, i.e., a fitted Energy Tree.
#' @param main Optional title for the plot.
#' @param type Character specifying the complexity of the plot:
#'   \code{extended} tries to visualize the distribution of the response
#'   variable in each terminal node whereas \code{simple} only gives some
#'   summary information.
#' @param terminal_panel Optional panel function of the form
#'   \code{function(node)} plotting the terminal nodes. Alternatively, a panel
#'   generating function of class \code{"grapcon_generator"} that is called with
#'   arguments \code{x} and \code{tp_args} to set up a panel function. By
#'   default, an appropriate panel function is chosen depending on the scale of
#'   the dependent variable.
#' @param tp_args	List of arguments passed to \code{terminal_panel} if this is
#'   a \code{"grapcon_generator"} object.
#' @param inner_panel	Optional panel function of the form
#'   \code{function(node)} plotting the inner nodes. Alternatively, a panel
#'   generating function of class \code{"grapcon_generator"} that is called with
#'   arguments \code{x} and \code{ip_args} to set up a panel function.
#' @param ip_args	List of arguments passed to \code{inner_panel} if this is a
#'   \code{"grapcon_generator"} object.
#' @param edge_panel Optional panel function of the form
#'   \code{function(split, ordered = FALSE, left = TRUE)} plotting the edges.
#'   Alternatively, a panel generating function of class
#'   \code{"grapcon_generator"} that is called with arguments \code{x} and
#'   \code{ep_args} to set up a panel function.
#' @param ep_args	List of arguments passed to \code{edge_panel} if this is a
#'   \code{"grapcon_generator"} object.
#' @param drop_terminal	Logical indicating whether all terminal nodes should
#'   be plotted at the bottom.
#' @param tnex Numeric value giving the terminal node extension in relation to
#'   the inner nodes.
#' @param newpage	Logical indicating whether \code{grid.newpage()} should be
#'   called.
#' @param pop	Logical indicating whether the viewport tree should be popped before
#'   return.
#' @param gp Graphical parameters.
#' @param ... Additional arguments.
#'
#' @details 
#' The \code{plot()} method for \code{"etree"} objects allows for the
#' visualization of fitted Energy Trees,  as returned by
#' \code{\link[etree:etree]{etree()}} or as contained in the \code{ensemble}
#' element of a fitted Random Energy Forest.
# The method is built upon the
# \code{partykit} version, from which it inherits several useful
# functionalities (cf. \code{\link[partykit]{party-plot}}).
#'
# @examples
#'
#' @returns 
#' No return value, called for side effects (plotting the tree).
#'
#' @method plot etree
#' @export

plot.etree <- function(x, main = NULL,
                       terminal_panel = NULL, tp_args = list(),
                       inner_panel = node_inner, ip_args = list(),
                       edge_panel = edge_simple, ep_args = list(),
                       type = c("extended", "simple"), drop_terminal = NULL, tnex = NULL,
                       newpage = TRUE, pop = TRUE, gp = gpar(), ...)
{
  ### compute default settings
  type <- match.arg(type)
  if (type == "simple") {
    x <- as.simpleparty(x)
    if (is.null(terminal_panel))
      terminal_panel <- node_terminal
    if (is.null(tnex)) tnex <- 1
    if (is.null(drop_terminal)) drop_terminal <- FALSE
    if (is.null(tp_args) || length(tp_args) < 1L) {
      tp_args <- list(FUN = .make_formatinfo_simpleparty(x, digits = getOption("digits") - 4L, sep = "\n"))
    } else {
      if(is.null(tp_args$FUN)) {
        tp_args$FUN <- .make_formatinfo_simpleparty(x, digits = getOption("digits") - 4L, sep = "\n")
      }
    }
  } else {
    if (is.null(terminal_panel)) {
      cl <- class(x$fitted[["(response)"]])
      if("factor" %in% cl) {
        terminal_panel <- node_barplot
      } else if("Surv" %in% cl) {
        terminal_panel <- node_surv
      } else if ("data.frame" %in% cl) {
        terminal_panel <- node_mvar
        if (is.null(tnex)) tnex <- 2 * NCOL(x$fitted[["(response)"]])
      } else {
        terminal_panel <- node_boxplot
      }
    }
    if (is.null(tnex)) tnex <- 2
    if (is.null(drop_terminal)) drop_terminal <- TRUE
  }
  
  plot.party(x, main = main,
             terminal_panel = terminal_panel, tp_args = tp_args,
             inner_panel = inner_panel, ip_args = ip_args,
             edge_panel = edge_panel, ep_args = ep_args,
             drop_terminal = drop_terminal, tnex = tnex,
             newpage = newpage, pop = pop, gp = gp, ...)
}
