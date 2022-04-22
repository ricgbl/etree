#' Predictions for Energy Trees
#' 
#' Compute predictions for objects of class \code{"etree"} (i.e., fitted Energy
#' Trees as returned by \code{\link[etree:etree]{etree()}}, or as contained in
#' the \code{ensemble} element of a fitted Random Energy Forest).
#'
#' @param object A fitted Energy Tree of class \code{"etree"}.
#' @param newdata Optional set of new covariates used to make predictions.
#'   Must be provided as a list, where each element is a different variable.
#'   Currently available types and the form they need to have to be correctly
#'   recognized are the following:
#' \itemize{
#'   \item Numeric: numeric or integer vectors;
#'   \item Nominal: factors;
#'   \item Functions: objects of class \code{"fdata"};
#'   \item Graphs: (lists of) objects of class \code{"igraph"}.
#   \item Persistence diagrams: (lists of) objects with
#   \code{attributes(x)$names == "diagram"}.
#' } 
#' Each element (i.e., variable) in the covariates list must have the same
#' \code{length()}, which corresponds to the (new) sample size. If
#' \code{newdata} is omitted, fitted values are returned.
#' @param perm Optional character vector of variable names. Splits of nodes
#'   with a primary split in any of these variables will be permuted (after
#'   dealing with surrogates). Note that surrogate split in the \code{perm}
#'   variables will not be permuted.
#' @param ... Additional arguments.
#'
#' @details
#' The \code{predict()} method for \code{"etree"} objects yields predictions for
#' fitted Energy Trees as returned by \code{\link[etree:etree]{etree()}} or as
#' contained in the \code{ensemble} element of a fitted Random Energy Forest.
#' Predictions are based either on fitted values (if \code{newdata} is
#' \code{NULL}) or on the new set of covariates (if \code{newdata} is provided).
#' The values of \code{split_type} and \code{coeff_split_type}, as well as the
#' number of components for each structured covariate (needed to compute an
#' equivalent representation for the covariates in \code{newdata} when
#' \code{split_type = "coeff"}), are automatically retrieved from the object of
#' class \code{"etree"}.
#' 
#' @section Value:
#' Predictions, in the form of a factor for classification or as a numeric
#' vector for regression.
#' 
# @examples
#'  
#' @method predict etree
#' @export

predict.etree <- function(object, newdata = NULL, perm = NULL, ...){
  
  # Retrieve split type
  split_type <- attr(object, 'split_type')
  
  # Retrieve terminal and inner nodes' ids, plus all basid in the 'coeff' case
  terminal <- nodeids(object, terminal = TRUE)
  if(max(terminal) > 1L) {
    inner <- 1L:max(terminal)
    inner <- inner[-terminal]
    # if(split_type == "coeff"){
    #   basids <<- nodeapply(object, ids = inner, by_node = TRUE,
    #                       FUN = function(node) basid_split(split_node(node)))
    #   train_max_basid <- max(unlist(basids))
    # }
  }
  
  # Coefficient expansion in the 'coeff' case
  if (!is.null(newdata) && split_type == "coeff") {
    
    cov_numbasis <- lapply(object$data,
                           function(cov) {
                             bases_names <- colnames(cov)
                             #colnames(), not ncol(), bc graph_shell may drop cols
                             if (!is.null(bases_names)) {
                               numbasis <- max(as.integer(bases_names))
                             }
                           })
    #number of basis for fda, max shell index for graphs, NULL for others
    #FIXME: only works if covs ordering in newdata is preserved (wrt to train)
    newdata <- mapply(function(cov, numbas) {
      attr(cov, 'numbasis') <- numbas
      return(cov)
    },
    newdata,
    cov_numbasis,
    SIMPLIFY = FALSE)
    
    newdata <- lapply(newdata, function(j) {
      
      switch(class(j),
             
             fdata = {
               
               train_nb <- attr(j, 'numbasis')
               fdata_est <- fda.usc::optim.basis(j,
                                                 numbasis = train_nb)
               coefs <- fda.usc::fdata2fd(fdata_est$fdata.est,
                                          type.basis = "bspline",
                                          nbasis = fdata_est$numbasis.opt)$coefs
               newcov <- data.frame(t(coefs))
               names(newcov) <- 1:length(names(newcov))
               return(newcov)
               
             },
             
             list = if (all(sapply(j, class) == 'igraph')) {
               
               train_nb <- attr(j, 'numbasis')
               newcov <- graph_shell(j,
                                     predicting = TRUE,
                                     max_shell = train_nb)
               return(newcov)
               
             },
             
             return(j)
             
      )})
  }
  
  ### compute fitted node ids first
  fitted <- if(is.null(newdata) && is.null(perm)) {
    object$fitted[["(fitted)"]]
  } else {
    if (is.null(newdata)) newdata <- model.frame(object)
    ### make sure all the elements in newdata have the same number of rows
    stopifnot(length(unique(sapply(newdata, NROW))) == 1L)
    
    if(max(terminal) == 1L) {
      rep.int(1L, unique(sapply(newdata, NROW)))
    } else {
      
      primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
        varid_split(split_node(node))
      })
      surrogate_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
        surr <- surrogates_node(node)
        if(is.null(surr)) return(NULL) else return(sapply(surr, varid_split))
      })
      vnames <- names(object$data)
      
      ### the splits of nodes with a primary split in perm
      ### will be permuted
      if (!is.null(perm)) {
        if (is.character(perm)) {
          stopifnot(all(perm %in% vnames))
          perm <- match(perm, vnames)
        } else {
          ### perm is a named list of factors coding strata
          ### (for varimp(..., conditional = TRUE)
          stopifnot(all(names(perm) %in% vnames))
          stopifnot(all(sapply(perm, is.factor)))
          tmp <- vector(mode = "list", length = length(vnames))
          tmp[match(names(perm), vnames)] <- perm
          perm <- tmp
        }
      }
      
      ## ## FIXME: the is.na() call takes loooong on large data sets
      ## unames <- if(any(sapply(newdata, is.na)))
      ##     vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
      ## else
      ##     vnames[unique(unlist(primary_vars))]
      unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
      
      vclass <- structure(lapply(object$data, class), .Names = vnames)
      ndnames <- names(newdata)
      ndclass <- structure(lapply(newdata, class), .Names = ndnames)
      checkclass <- all(sapply(unames, function(x)
        isTRUE(all.equal(vclass[[x]], ndclass[[x]]))))
      factors <- sapply(unames, function(x) inherits(object$data[[x]], "factor"))
      checkfactors <- all(sapply(unames[factors], function(x)
        isTRUE(all.equal(levels(object$data[[x]]), levels(newdata[[x]])))))
      ## FIXME: inform about wrong classes / factor levels?
      if(all(unames %in% ndnames) && checkclass && checkfactors) {
        vmatch <- match(vnames, ndnames)
        fitted_node_predict(node_party(object), data = newdata,
                            vmatch = vmatch, perm = perm)
      } else {
        if (!is.null(object$terms)) {
          ### <FIXME> this won't work for multivariate responses
          ### </FIXME>
          xlev <- lapply(unames[factors],
                         function(x) levels(object$data[[x]]))
          names(xlev) <- unames[factors]
          #         mf <- model.frame(delete.response(object$terms), newdata,
          #                          xlev = xlev)
          # fitted_node_predict(node_party(object), data = newdata,
          #             vmatch = match(vnames, names(mf)), perm = perm)
          fitted_node_predict(node_party(object), data = newdata,
                              perm = perm)
        } else
          stop("") ## FIXME: write error message
      }
    }
  }
  ### compute predictions
  predict_party(object, fitted, newdata, ...)
}
