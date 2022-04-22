#' Methods for "etree" objects
#' 
#' Methods for objects of class \code{"etree"}.
#' 
#' @aliases print.etree length.etree depth.etree width.etree "[.etree" "[[.etree"
#' 
#' @param x Object of class \code{"etree"}.
#' @param i Integer specifying the root of the subtree to extract.
# @param terminal_panel Panel function for printing terminal nodes.
# @param tp_args List containing arguments to \code{terminal_panel}.
# @param inner_panel Panel function for printing inner nodes.
# @param ip_args List containing argument to \code{inner_panel}.
# @param header_panel Panel function for printing the header.
# @param footer_panel Panel function for printing the footer.
#' @param FUN Function to be applied to nodes.
#' @param digits Number of digits to be printed.
#' @param header Header to be printed.
#' @param footer Footer to be printed.
#' @param root Logical indicating whether the root node should be counted in
#'   \code{depth()} or not (default).
#' @param ... Additional arguments.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' ## Covariates
#' nobs <- 100
#' cov_num <- rnorm(nobs)
#' cov_nom <- factor(rbinom(nobs, size = 1, prob = 0.5))
#' cov_gph <- lapply(1:nobs, function(j) igraph::sample_gnp(100, 0.2))
#' cov_fun <- fda.usc::rproc2fdata(nobs, seq(0, 1, len = 100), sigma = 1)
# x <- lapply(rep(100, nobs), function(np) TDA::circleUnif(np))
# cov_per <- lapply(x, TDA::ripsDiag, maxdimension = 1, maxscale = 3)  
#' cov_list <- list(cov_num, cov_nom, cov_gph, cov_fun)
#' 
#' ## Response variable
#' resp_reg <- cov_num ^ 2
#' 
#' ## Fit
#' etree_fit <- etree(response = resp_reg, covariates = cov_list)
#' 
#' ## Print
#' print(etree_fit)
#' 
#' ## Number of nodes in the tree
#' length(etree_fit)
#' 
#' ## Depth of the tree
#' depth(etree_fit)
#' 
#' ## Number of terminal nodes in the tree
#' width(etree_fit)
#' 
#' ## Extract subtrees
#' etree_fit[2]
#' etree_fit[[2]]
#' 
#' }
#' 
#' @name etree-methods
NULL


#' @describeIn etree-methods Generates textual representation of the tree.
#' @method print etree
#' @export
print.etree <- function(x,
                        FUN = NULL, digits = getOption("digits") - 4,
                        header = NULL, footer = TRUE, ...)
{
  if(is.null(FUN)) return(print(as.simpleparty(x), digits = digits,
                                header = header, footer = footer, ...))
  
  digits <- max(c(0, digits))
  
  ## FIXME: terms/call/? for "ctree" objects
  if(is.null(header)) header <- !is.null(terms(x))
  header_panel <- if(header) function(party) {
    c("", "Model formula:", deparse(formula(terms(party))), "", "Fitted party:", "")
  } else function(party) ""
  
  footer_panel <- if(footer) function(party) {
    n <- width(party)
    n <- format(c(length(party) - n, n))
    
    c("", paste("Number of inner nodes:   ", n[1]),
      paste("Number of terminal nodes:", n[2]), "")
  } else function (party) ""
  
  y <- x$fitted[["(response)"]]
  w <- x$fitted[["(weights)"]]
  if(is.null(w)) {
    wdigits <- 0
    wsym <- "n"
  } else {
    if(isTRUE(all.equal(w, round(w)))) {
      wdigits <- 0
      wsym <- "n"
    } else {
      wdigits <- max(c(0, digits - 2))
      wsym <- "w"
    }
  }
  yclass <- class(y)[1]
  if(yclass == "ordered") yclass <- "factor"
  if(!(yclass %in% c("Surv", "factor"))) yclass <- "numeric"
  
  if(is.null(FUN)) FUN <- switch(yclass,
                                 "numeric" = function(y, w, digits) {
                                   yhat <- .pred_numeric_response(y, w)
                                   yerr <- sum(w * (y - yhat)^2)
                                   digits2 <- max(c(0, digits - 2))
                                   paste(format(round(yhat, digits = digits), nsmall = digits),
                                         " (", wsym, " = ", format(round(sum(w), digits = wdigits), nsmall = wdigits), ", err = ",
                                         format(round(yerr, digits = digits2), nsmall = digits2), ")", sep = "")
                                 },
                                 "Surv" = function(y, w, digits) {
                                   paste(format(round(.pred_Surv_response(y, w), digits = digits), nsmall = digits),
                                         " (", wsym, " = ", format(round(sum(w), digits = wdigits), nsmall = wdigits), ")", sep = "")
                                 },
                                 "factor" = function(y, w, digits) {
                                   tab <- round(.pred_factor(y, w) * sum(w))
                                   mc <- round(100 * (1 - max(tab)/sum(w)), digits = max(c(0, digits - 2)))
                                   paste(names(tab)[which.max(tab)], " (",
                                         wsym, " = ", format(round(sum(w), digits = wdigits), nsmall = wdigits),
                                         ", err = ", mc, "%)", sep = "")
                                 }
  )
  
  node_labs <- nodeapply(x, nodeids(x), function(node) {
    y1 <- node$fitted[["(response)"]]
    w <- node$fitted[["(weights)"]]
    if(is.null(w)) w <- rep.int(1L, NROW(y1))
    FUN(y1, w, digits)
  }, by_node = FALSE)
  node_labs <- paste(":", format(do.call("c", node_labs)))
  
  terminal_panel <- function(node) node_labs[id_node(node)]
  
  print.party(x, terminal_panel = terminal_panel,
              header_panel = header_panel, footer_panel = footer_panel, ...)
}


#' @describeIn etree-methods Number of nodes in the tree.
#' @method length etree
#' @export
length.party <- length.etree <- function(x) {
  length(nodeids(x))
}


#' @describeIn etree-methods Depth of the three.
#' @method depth etree
#' @export
depth.party <- depth.etree <- function(x, root = FALSE, ...) {
  depth(node_party(x), root = root, ...)
}


#' @describeIn etree-methods Number of terminal nodes.
#' @method width etree
#' @export
width.party <- width.etree <- function(x, ...) {
  width(node_party(x), ...)
}


#' @describeIn etree-methods Extract subtrees.
#' @usage \method{[}{etree}(x, i, \dots)
#' @method "[" etree
#' @export
"[.etree" <- function(x, i, ...) {
  if (is.character(i) && !is.null(names(x)))
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[,findx:ncol(dat), drop = FALSE]
    dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0)
      dat <- x$data
  } else {
    fit <- NULL
    dat <- x$data
  }
  nam <- names(x)[nodeids(x, from = i, terminal = FALSE)]
  
  recFun <- function(node) {
    if (id_node(node) == i) return(node)
    kid <- sapply(kids_node(node), id_node)
    return(recFun(node[[max(which(kid <= i))]]))
  }
  node <- recFun(node_party(x))
  names(dat) <- names(x$data)
  #kind of hack, but needed to print the correct variables when subsetting
  #TODO: find a better way to solve this
  ret <- party(node = node, data = dat, fitted = fit, 
               terms = x$terms, names = nam, info = x$info)
  class(ret) <- class(x)
  ret
}

#' @describeIn etree-methods Extract subtrees.
#' @usage \method{[[}{etree}(x, i, \dots)
#' @method "[[" etree
#' @export
"[[.etree" <- function(x, i, ...) {
  if (is.character(i) && !is.null(names(x)))
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[,findx:ncol(dat), drop = FALSE]
    dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0)
      dat <- x$data
  } else {
    fit <- NULL
    dat <- x$data
  }
  nam <- names(x)[nodeids(x, from = i, terminal = FALSE)]
  
  recFun <- function(node) {
    if (id_node(node) == i) return(node)
    kid <- sapply(kids_node(node), id_node)
    return(recFun(node[[max(which(kid <= i))]]))
  }
  node <- recFun(node_party(x))
  names(dat) <- names(x$data)
  #kind of hack, but needed to print the correct variables when subsetting
  #TODO: find a better way to solve this
  ret <- party(node = node, data = dat, fitted = fit, 
               terms = x$terms, names = nam, info = x$info)
  class(ret) <- class(x)
  ret
}

