#' Apply functions over nodes
#' 
#' Returns a list of values obtained by applying a function to \code{"etree"} or
#' \code{"partynode"} objects.
#' 
#' @aliases nodeapply.partynode nodeapply.etree
#'
#' @param obj Object of class \code{"etree"} or \code{"partynode"}.
#' @param ids Integer vector of node identifiers to apply over.
#' @param FUN Function to be applied to nodes. By default, the node itself is
#'   returned.
#' @param by_node Logical indicating whether FUN should be applied to subsets of
#'   \code{"partynode"} objects (default) or not, in which case it is applied to
#'   subsets of \code{"etree"} objects.
#' @param ... Additional arguments.
#' 
#' @section Value:
#' A list of results whose length is given by \code{length(ids)}.
#' 
#' @details 
#' The method for \code{"partynode"} objects apply function \code{FUN} to all
#' nodes with node identifiers in \code{ids}. The method for \code{"etree"}
#' objects by default calls the \code{nodeapply} method on the corresponding
#' node slot. If \code{by_node} is \code{FALSE}, it is applied to the
#' \code{"etree"} object with root node \code{ids}.
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
#' ## Get pvalues of inner nodes
#' tnodes <- nodeids(etree_fit, terminal = TRUE)
#' nodes <- 1:max(tnodes)
#' inodes <- nodes[-tnodes]
#' nodeapply(etree_fit, ids = inodes, FUN = function(n) n$info$pvalue)
#' }

#' @export
nodeapply <- function(obj, ids = 1, FUN = NULL, ...)
  UseMethod("nodeapply")


#' @describeIn nodeapply nodeapply() method for objects of class "partynode".
#' @method nodeapply partynode
nodeapply.partynode <- function(obj, ids = 1, FUN = NULL, ...) {
  
  stopifnot(isTRUE(all.equal(ids, round(ids))))
  ids <- as.integer(ids)
  
  if(is.null(FUN)) FUN <- function(x, ...) x
  
  if (length(ids) == 0)
    return(NULL)
  
  rval <- vector(mode = "list", length = length(ids))
  rval_id <- rep(0, length(ids))
  i <- 1
  
  recFUN <- function(node, ...) {
    if(id_node(node) %in% ids) {
      rval_id[i] <<- id_node(node)
      rval[[i]] <<- FUN(node, ...)
      i <<- i + 1
    }
    kids <- kids_node(node)
    if(length(kids) > 0) {
      for(j in 1:length(kids)) recFUN(kids[[j]])
    }
    invisible(TRUE)
  }
  foo <- recFUN(obj)
  rval <- rval[match(ids, rval_id)]
  return(rval)
}

#' @describeIn nodeapply nodeapply() method for objects of class "etree".
#' @method nodeapply etree
#' @export
nodeapply.party <- nodeapply.etree <-
  function(obj, ids = 1, FUN = NULL, by_node = TRUE, ...) {
    
    stopifnot(isTRUE(all.equal(ids, round(ids))))
    ids <- as.integer(ids)
    
    if(is.null(FUN)) FUN <- function(x, ...) x
    
    if (length(ids) == 0)
      return(NULL)
    
    if (by_node) {
      rval <- nodeapply(node_party(obj), ids = ids, FUN = FUN, ...)
    } else {
      rval <- lapply(ids, function(i) FUN(obj[[i]], ...))
    }
    
    names(rval) <- names(obj)[ids]
    return(rval)
  }



#' Extract node identifiers.
#'
#' Extract unique identifiers from inner and terminals nodes of \code{"etree"}
#' or \code{"partynode"} objects.
#'
#' @aliases nodeids.partynode nodeids.etree
#'
#' @param obj Object of class \code{"etree"} or \code{"partynode"}.
#' @param from Integer specifying the node to start from.
#' @param terminal Logical indicating whether only identifiers of terminal nodes
#'   should be returned (\code{FALSE} by default).
#' @param ... Additional arguments.
#' 
#' @section Value:
#' An integer vector of node identifiers.
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
#' ## Get all nodes identifiers
#' nodes_ids <- nodeids(etree_fit)
#' 
#' ## Get terminal nodes identifiers
#' tnodes_ids <- nodeids(etree_fit, terminal = TRUE)
#' 
#' ## Get all nodes identifiers starting from 2
#' nodes_ids2 <- nodeids(etree_fit, from = 2)
#' }

#' @export
nodeids <- function(obj, ...)
  UseMethod("nodeids")

#' @describeIn nodeids nodeids() method for objects of class "partynode".
#' @method nodeids partynode
nodeids.partynode <- function(obj, from = NULL, terminal = FALSE, ...) {
  
  if(is.null(from)) from <- id_node(obj)
  
  id <- function(node, record = TRUE, terminal = FALSE) {
    if(!record) return(NULL)
    if(!terminal)
      return(id_node(node))
    else
      if(is.terminal(node)) return(id_node(node)) else return(NULL)
  }
  
  rid <- function(node, record = TRUE, terminal = FALSE) {  
    myid <- id(node, record = record, terminal = terminal)
    if(is.terminal(node)) return(myid)
    kids <- kids_node(node)    
    kids_record <- if(record)  
      rep(TRUE, length(kids))
    else
      sapply(kids, id_node) == from
    return(c(myid,
             unlist(lapply(1:length(kids), function(i)
               rid(kids[[i]], record = kids_record[i], terminal = terminal)))
    ))
  }
  
  return(rid(obj, from == id_node(obj), terminal))
}

#' @describeIn nodeids nodeids() method for objects of class "etree".
#' @method nodeids etree
#' @export
nodeids.party <- nodeids.etree <-
  function(obj, from = NULL, terminal = FALSE, ...)
    nodeids(node_party(obj), from = from, terminal = terminal, ...)

