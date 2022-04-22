#' etree: Classification and Regression With Structured and Mixed-Type Data
#' 
#' Implementation of Energy Trees, a statistical model to perform classification
#' and regression with structured and mixed-type data. This is made possible by
#' starting from the statistically sound structure of Conditional Trees, and
#' using Energy Statistics to test independence between data objects of possibly
#' structured and possibly different nature. The type of structured objects
#' covered so far are functions and graphs. The package builds upon 'partykit'
#' to provide functionalities for fitting, printing, plotting, and predicting
#' with Energy Trees.
#' 
#' @rawNamespace import(partykit, except = c(partynode, kidids_node,
#'   fitted_node, party, predict, predict_party, edge_simple, `[`, `[[`,
#'   partysplit, kidids_split, node_barplot, node_boxplot, node_surv, node_ecdf,
#'   node_mvar, data_party, node_inner, nodeapply, nodeids, print, depth,
#'   width))
#' @importFrom grDevices gray.colors
#' @importFrom graphics boxplot
#' @importFrom stats approxfun as.dist density dist ecdf formula knots
#'   model.frame na.exclude p.adjust predict quantile terms var weighted.mean
#' @importFrom utils combn
#' @importFrom grid grid.clip grid.rect grid.layout grid.lines grid.points
#' grid.polygon grid.text grid.xaxis grid.yaxis gpar depth popViewport
#' pushViewport unit upViewport viewport
#' @importFrom survival survfit
#'
#' @docType package
#' @name etree-package

NULL
