
# node --------------------------------------------------------------------

partynode <- function(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL, centroids = NULL) {
  
  if (!is.integer(id) || length(id) != 1) {
    id <- as.integer(id0 <- id)
    if (any(is.na(id)) || !isTRUE(all.equal(id0, id)) || length(id) != 1)
      stop(sQuote("id"), " ", "must be a single integer")
  }
  
  if (is.null(split) != is.null(kids)) {
    stop(sQuote("split"), " ", "and", " ", sQuote("kids"), " ",
         "must either both be specified or unspecified")
  }
  
  if (!is.null(kids)) {
    if (!(is.list(kids) && all(sapply(kids, inherits, "partynode")))
        || length(kids) < 2)
      stop(sQuote("kids"), " ", "must be an integer vector or a list of",
           " ", sQuote("partynode"), " ", "objects")
  }
  
  if (!is.null(surrogates)) {
    if (!is.list(surrogates) || any(!sapply(surrogates, inherits, "partysplit")))
      stop(sQuote("split"), " ", "is not a list of", " ", sQuote("partysplit"),
           " ", "objects")
  }
  
  node <- list(id = id, split = split, kids = kids, surrogates = surrogates, info = info, centroids = centroids)
  class(node) <- "partynode"
  return(node)
}

as.partynode.list <- function(x, ...) {
  
  if (!all(sapply(x, inherits, what = "list")))
    stop("'x' has to be a list of lists")
  
  if (!all(sapply(x, function(x) "id" %in% names(x))))
    stop("each list in 'x' has to define a node 'id'")
  
  ok <- sapply(x, function(x)
    all(names(x) %in% c("id", "split", "kids", "surrogates", "info", "centroids")))
  if (any(!ok))
    sapply(which(!ok), function(i)
      warning(paste("list element", i, "defines additional elements:",
                    paste(names(x[[i]])[!(names(x[[i]]) %in%
                                            c("id", "split", "kids", "surrogates", "info", "centroids"))],
                          collapse = ", "))))
  
  ids <- as.integer(sapply(x, function(node) node$id))
  if(any(duplicated(ids))) stop("nodeids must be unique integers")
  x <- x[order(ids)]
  ids <- ids[order(ids)]
  
  new_recnode <- function(i) {
    x_i <- x[[which(ids == i)]]
    if (is.null(x_i$kids))
      partynode(id = x_i$id, info = x_i$info)
    else
      partynode(id = x_i$id, split = x_i$split,
                kids = lapply(x_i$kids, new_recnode),
                surrogates = x_i$surrogates,
                info = x_i$info)
  }
  
  ret <- new_recnode(ids[1L])
  ### <FIXME> duplicates recursion but makes sure
  ###    that the ids are in pre-order notation with
  ###    from defined in as.partynode.partynode
  ### </FIXME>
  as.partynode(ret, ...)
}

kidids_node <- function(node, data, vmatch = 1:length(data), obs = NULL,
                        perm = NULL) {
  
  primary <- split_node(node)
  surrogates <- surrogates_node(node)
  
  ### perform primary split
  x <- kidids_split(primary, data, vmatch, obs)
  
  ### surrogate / random splits if needed
  if (any(is.na(x))) {
    ### surrogate splits
    if (length(surrogates) >= 1) {
      for (surr in surrogates) {
        nax <- is.na(x)
        if (!any(nax)) break;
        x[nax] <- kidids_split(surr, data, vmatch, obs = obs)[nax]
      }
    }
    nax <- is.na(x)
    ### random splits
    if (any(nax)) {
      prob <- prob_split(primary)
      x[nax] <- sample(1:length(prob), sum(nax), prob = prob,
                       replace = TRUE)
    }
  }
  
  ### permute variable `perm' _after_ dealing with surrogates etc.
  if (!is.null(perm)) {
    if (is.integer(perm)) {
      if (varid_split(primary) %in% perm)
        x <- .resample(x)
    } else {
      if (is.null(obs)) obs <- 1:length(data)
      strata <- perm[[varid_split(primary)]]
      if (!is.null(strata)) {
        strata <- strata[obs, drop = TRUE]
        for (s in levels(strata))
          x[strata == s] <- .resample(x[strata == s])
      }
    }
  }
  return(x)
}

kidids_node_predict <- function(node, data, vmatch = 1:length(data), obs = NULL,
                                perm = NULL) {
  
  primary <- split_node(node)
  surrogates <- surrogates_node(node)
  
  ### perform primary split
  x <- kidids_split_predict(primary, data, vmatch, obs)
  
  ### surrogate / random splits if needed
  if (any(is.na(x))) {
    ### surrogate splits
    if (length(surrogates) >= 1) {
      for (surr in surrogates) {
        nax <- is.na(x)
        if (!any(nax)) break;
        x[nax] <- kidids_split_predict(surr, data, vmatch, obs = obs)[nax]
      }
    }
    nax <- is.na(x)
    ### random splits
    if (any(nax)) {
      prob <- prob_split(primary)
      x[nax] <- sample(1:length(prob), sum(nax), prob = prob,
                       replace = TRUE)
    }
  }
  
  ### permute variable `perm' _after_ dealing with surrogates etc.
  if (!is.null(perm)) {
    if (is.integer(perm)) {
      if (varid_split(primary) %in% perm)
        x <- .resample(x)
    } else {
      if (is.null(obs)) obs <- 1:length(data)
      strata <- perm[[varid_split(primary)]]
      if (!is.null(strata)) {
        strata <- strata[obs, drop = TRUE]
        for (s in levels(strata))
          x[strata == s] <- .resample(x[strata == s])
      }
    }
  }
  return(x)
}

fitted_node <- function(node, data, vmatch = 1:length(data),
                        obs = 1:unique(sapply(data, NROW)), perm = NULL) {
  
  if (is.logical(obs)) obs <- which(obs)
  if (is.terminal(node))
    return(rep(id_node(node), length(obs)))
  retid <- nextid <- kidids_node(node, data, vmatch, obs, perm)
  
  for (i in unique(nextid)) {
    indx <- nextid == i
    retid[indx] <- fitted_node(kids_node(node)[[i]], data,
                               vmatch, obs[indx], perm)
  }
  return(retid)
}

fitted_node_predict <- function(node, data, vmatch = 1:length(data),
                                obs = 1:unique(sapply(data, NROW)), perm = NULL) {
  
  if (is.logical(obs)) obs <- which(obs)
  if (is.terminal(node))
    return(rep(id_node(node), length(obs)))
  retid <- nextid <- kidids_node_predict(node, data, vmatch, obs, perm)
  
  for (i in unique(nextid)) {
    indx <- nextid == i
    retid[indx] <- fitted_node_predict(kids_node(node)[[i]], data,
                                       vmatch, obs[indx], perm)
  }
  return(retid)
}

is.terminal.partynode <- function(x, ...) {
  kids <- is.null(kids_node(x))
  split <- is.null(split_node(x))
  if (kids != split)
    stop("x", " ", "is incorrect node")
  kids
}

#' @method depth partynode
#' @export
depth.partynode <- function(x, root = FALSE, ...) {
  if (is.terminal(x)) return(as.integer(root))
  max(sapply(kids_node(x), depth, root = root)) + 1L
}

width.partynode <- function(x, ...) {
  if (is.terminal(x)) return(1)
  sum(sapply(kids_node(x), width.partynode))
}


# party -------------------------------------------------------------------

party <- function(node, data, fitted = NULL, terms = NULL, names = NULL, info = NULL) {
  
  stopifnot(inherits(node, "partynode"))
  stopifnot(inherits(data, "list") | inherits(data, "data.frame"))
  #without adding 'data.frame', it throws error for CLS plots
  ### make sure all split variables are there
  ids <- nodeids(node)[!nodeids(node) %in% nodeids(node, terminal = TRUE)]
  varids <- unique(unlist(nodeapply(node, ids = ids, FUN = function(x)
    varid_split(split_node(x)))))
  stopifnot(varids %in% 1:length(data))

  if(!is.null(fitted)) {
    stopifnot(inherits(fitted, "data.frame"))
    stopifnot(all(sapply(data, NROW) == 0L) | all(sapply(data, NROW) == NROW(fitted)))
    
    # try to provide default variable "(fitted)"
    if(all(sapply(data, NROW) > 0L)) {
      if(!("(fitted)" %in% names(fitted)))
        fitted[["(fitted)"]] <- fitted_node(node, data = data)
    } else {
      stopifnot("(fitted)" == names(fitted)[1L])
    }

    nt <- nodeids(node, terminal = TRUE)
    stopifnot(all(fitted[["(fitted)"]] %in% nt))
    
    node <- as.partynode(node, from = 1L)
    nt2 <- nodeids(node, terminal = TRUE)
    fitted[["(fitted)"]] <- nt2[match(fitted[["(fitted)"]], nt)]
  } else {
    node <- as.partynode(node, from = 1L)
    # default "(fitted)"
    if(all(sapply(data, NROW) > 0L) & missing(fitted))
      fitted <- data.frame("(fitted)" = fitted_node(node,
                                                    data = data), check.names = FALSE)
  }

  party <- list(node = node, data = data, fitted = fitted,
                terms = NULL, names = NULL, info = info)
  class(party) <- "party"

  if(!is.null(terms)) {
    stopifnot(inherits(terms, "terms"))
    party$terms <- terms
  }
  
  if (!is.null(names)) {
    n <- length(nodeids(party, terminal = FALSE))
    if (length(names) != n)
      stop("invalid", " ", sQuote("names"), " ", "argument")
    party$names <- names
  }
  
  party
}

.names_party <- function(party) {
  names <- party$names
  if (is.null(names))
    names <- as.character(nodeids(party, terminal = FALSE))
  names
}

predict_party <- function(party, id, newdata = NULL, ...)
  UseMethod("predict_party")

### do nothing except returning the fitted ids
predict_party.default <- function(party, id, newdata = NULL, FUN = NULL, ...) {
  
  if (length(list(...)) > 1)
    warning("argument(s)", " ", sQuote(names(list(...))), " ", "have been ignored")
  
  ## get observation names: either node names or
  ## observation names from newdata
  nam <- if(is.null(newdata)) {
    if(is.null(rnam <- rownames(data_party(party)))) names(party)[id] else rnam
  } else {
    rownames(newdata[[1]])
  }
  if(length(nam) != length(id)) nam <- NULL
  
  if (!is.null(FUN))
    return(.simplify_pred(nodeapply(party,
                                    nodeids(party, terminal = TRUE), FUN, by_node = TRUE), id, nam))
  
  ## special case: fitted ids
  return(structure(id, .Names = nam))
}

predict_party.constparty <- function(party, id, newdata = NULL,
                                     type = c("response", "prob", "quantile", "density", "node"),
                                     at = if (type == "quantile") c(0.1, 0.5, 0.9),
                                     FUN = NULL, simplify = TRUE, ...)
{
  ## extract fitted information
  response <- party$fitted[["(response)"]]
  weights <- party$fitted[["(weights)"]]
  fitted <- party$fitted[["(fitted)"]]
  if (is.null(weights)) weights <- rep(1, NROW(response))
  
  ## get observation names: either node names or
  ## observation names from newdata
  nam <- if(is.null(newdata)) names(party)[id] else rownames(newdata[[1]])
  if(length(nam) != length(id)) nam <- NULL
  
  ## match type
  type <- match.arg(type)
  
  ## special case: fitted ids
  if(type == "node")
    return(structure(id, .Names = nam))
  
  ### multivariate response
  if (is.data.frame(response)) {
    ret <- lapply(response, function(r) {
      ret <- .predict_party_constparty(node_party(party), fitted = fitted,
                                       response = r, weights, id = id, type = type, at = at, FUN = FUN, ...)
      if (simplify) .simplify_pred(ret, id, nam) else ret
    })
    if (all(sapply(ret, is.atomic)))
      ret <- as.data.frame(ret)
    names(ret) <- colnames(response)
    return(ret)
  }
  
  ### univariate response
  ret <- .predict_party_constparty(node_party(party), fitted = fitted, response = response,
                                   weights = weights, id = id, type = type, at = at, FUN = FUN, ...)
  if (simplify) .simplify_pred(ret, id, nam) else ret[as.character(id)]
}

### functions for node prediction based on fitted / response
.pred_Surv <- function(y, w) {
  if (length(y) == 0) return(NA)
  survfit(y ~ 1, weights = w, subset = w > 0)
}

.pred_Surv_response <- function(y, w) {
  if (length(y) == 0) return(NA)
  .median_survival_time(.pred_Surv(y, w))
}

.pred_factor <- function(y, w) {
  lev <- levels(y)
  sumw <- tapply(w, y, sum)
  sumw[is.na(sumw)] <- 0
  prob <- sumw / sum(w)
  names(prob) <- lev
  return(prob)
}

.pred_factor_response <- function(y, w) {
  prob <- .pred_factor(y, w)
  return(factor(which.max(prob), levels = 1:nlevels(y),
                labels = levels(y),
                ordered = is.ordered(y)))
  return(prob)
}

.pred_numeric_response <- function(y, w)
  weighted.mean(y, w, na.rm = TRUE)

.pred_ecdf <- function(y, w) {
  if (length(y) == 0) return(NA)
  iw <- as.integer(round(w))
  if (max(abs(w - iw)) < sqrt(.Machine$double.eps)) {
    y <- rep(y, w)
    return(ecdf(y))
  } else {
    stop("cannot compute empirical distribution function with non-integer weights")
  }
}

.pred_quantile <- function(y, w) {
  y <- rep(y, w)
  function(p, ...) quantile(y, probs = p, ...)
}

.pred_density <- function(y, w) {
  d <- density(y, weights = w / sum(w))
  approxfun(d[1:2], rule = 2)
}

### workhorse: compute predictions based on fitted / response data
.predict_party_constparty <- function(node, fitted, response, weights,
                                      id = id, type = c("response", "prob", "quantile", "density"),
                                      at = if (type == "quantile") c(0.1, 0.5, 0.9), FUN = NULL, ...) {
  
  type <- match.arg(type)
  if (is.null(FUN)) {
    
    rtype <- class(response)[1]
    if (rtype == "ordered") rtype <- "factor"
    if (rtype == "integer") rtype <- "numeric"
    if (rtype == "AsIs") rtype <- "numeric"
    
    if (type %in% c("quantile", "density") && rtype != "numeric")
      stop("quantile and density estimation currently only implemented for numeric responses")
    
    FUN <- switch(rtype,
                  "Surv" = if (type == "response") .pred_Surv_response else .pred_Surv,
                  "factor" = if (type == "response") .pred_factor_response else .pred_factor,
                  "numeric" = switch(type,
                                     "response" = .pred_numeric_response,
                                     "prob" = .pred_ecdf,
                                     "quantile" = .pred_quantile,
                                     "density" = .pred_density)
    )
  }
  
  ## empirical distribution in each leaf
  if (all(id %in% fitted)) {
    tab <- tapply(1:NROW(response), fitted,
                  function(i) FUN(response[i], weights[i]), simplify = FALSE)
  } else {
    ### id may also refer to inner nodes
    tab <- as.array(lapply(sort(unique(id)), function(i) {
      index <- fitted %in% nodeids(node, i, terminal = TRUE)
      ret <- FUN(response[index], weights[index])
      ### no information about i in fitted
      if (all(!index)) ret[1] <- NA
      return(ret)
    }))
    names(tab) <- as.character(sort(unique(id)))
  }
  if (inherits(tab[[1]], "function") && !is.null(at))
    tab <- lapply(tab, function(f) f(at))
  tn <- names(tab)
  dim(tab) <- NULL
  names(tab) <- tn
  
  tab
}


### simplify structure of predictions
.simplify_pred <- function(tab, id, nam) {
  
  if (all(sapply(tab, length) == 1) & all(sapply(tab, is.atomic))) {
    ret <- do.call("c", tab)
    names(ret) <- names(tab)
    #ret <- if (is.factor(tab[[1]]))
    #  factor(ret[as.character(id)], levels = 1:length(levels(tab[[1]])),
    #        labels = levels(tab[[1]]), ordered = is.ordered(tab[[1]]))
    #else
    #  ret[as.character(id)]
    #names(ret) <- nam
    ret <- ret[as.character(id)]
  } else if (length(unique(sapply(tab, length))) == 1 &
             all(sapply(tab, is.numeric))) {
    ret <- matrix(unlist(tab), nrow = length(tab), byrow = TRUE)
    colnames(ret) <- names(tab[[1]])
    rownames(ret) <- names(tab)
    ret <- ret[as.character(id),, drop = FALSE]
    rownames(ret) <- nam
  } else {
    ret <- tab[as.character(id)]
    names(ret) <- nam
  }
  ret
}

data_party <- function(party, id = 1L)
  UseMethod("data_party")

data_party.default <- function(party, id = 1L) {
  
  extract <- function(id) {
    if(is.null(party$fitted))
      if(length(party$data) == 0) return(NULL)
    else
      stop("cannot subset data without fitted ids")
    
    ### which terminal nodes follow node number id?
    nt <- nodeids(party, id, terminal = TRUE)
    wi <- party$fitted[["(fitted)"]] %in% nt
    
    ret <- if (length(party$data) == 0)
      subset(party$fitted, wi)
    else
      subset(cbind(party$data, party$fitted), wi)
    ret
  }
  if (length(id) > 1)
    return(lapply(id, extract))
  else
    return(extract(id))
}

.list.rules.party <- function(x, i = NULL, ...) {
  if (is.null(i)) i <- nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    ret <- sapply(i, .list.rules.party, x = x)
    names(ret) <- if (is.character(i)) i else names(x)[i]
    return(ret)
  }
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
  
  rule <- c()
  
  recFun <- function(node) {
    if (id_node(node) == i) return(NULL)
    kid <- sapply(kids_node(node), id_node)
    whichkid <- max(which(kid <= i))
    split <- split_node(node)
    ivar <- varid_split(split)
    svar <- names(dat)[ivar]
    index <- index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index))
        index <- ((1:nlevels(dat[, svar])) > breaks_split(split)) + 1
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"",
                     paste(slevels, collapse = "\", \"", sep = ""), "\")",
                     sep = "")
    } else {
      if (is.null(index)) index <- 1:length(kid)
      breaks <- cbind(c(-Inf, breaks_split(split)),
                      c(breaks_split(split), Inf))
      sbreak <- breaks[index == whichkid,]
      right <- right_split(split)
      srule <- c()
      if (is.finite(sbreak[1]))
        srule <- c(srule,
                   paste(svar, ifelse(right, ">", ">="), sbreak[1]))
      if (is.finite(sbreak[2]))
        srule <- c(srule,
                   paste(svar, ifelse(right, "<=", "<"), sbreak[2]))
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  node <- recFun(node_party(x))
  paste(rule, collapse = " & ")
}


# plot --------------------------------------------------------------------

node_inner <- function(obj, id = TRUE, pval = TRUE, abbreviate = FALSE, fill = "white", gp = gpar())
{
  meta <- obj$data
  nam <- names(obj)
  
  extract_label <- function(node) {
    if(is.terminal(node)) return(rep.int("", 2L))
    
    varlab <- character_split(split_node(node), meta)$name
    varbas <- character_split(split_node(node), meta)$basid
    if(!is.null(varbas)) varlab <- paste0(varlab, '.', varbas)
    if(abbreviate > 0L) varlab <- abbreviate(varlab, as.integer(abbreviate))
    
    ## FIXME: make more flexible rather than special-casing p-value
    if(pval) {
      nullna <- function(x) is.null(x) || is.na(x)
      pval <- suppressWarnings(try(!nullna(info_node(node)$pvalue), silent = TRUE))
      pval <- if(inherits(pval, "try-error")) FALSE else pval
    }
    if(pval) {
      pvalue <- node$info$pvalue
      plab <- ifelse(pvalue < 10^(-3L),
                     paste("p <", 10^(-3L)),
                     paste("p =", round(pvalue, digits = 3L)))
    } else {
      plab <- ""
    }
    return(c(varlab, plab))
  }
  
  maxstr <- function(node) {
    lab <- extract_label(node)
    klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
    lab <- c(lab, klab)
    lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
    lab <- lab[which.max(nchar(lab))]
    if(length(lab) < 1L) lab <- ""
    return(lab)
  }
  
  nstr <- maxstr(node_party(obj))
  if(nchar(nstr) < 6) nstr <- "aAAAAa"
  
  ### panel function for the inner nodes
  rval <- function(node) {
    pushViewport(viewport(gp = gp, name = paste("node_inner", id_node(node), "_gpar", sep = "")))
    node_vp <- viewport(
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "strwidth", nstr) * 1.3,
      height = unit(3, "lines"),
      name = paste("node_inner", id_node(node), sep = ""),
      gp = gp
    )
    pushViewport(node_vp)
    
    xell <- c(seq(0, 0.2, by = 0.01),
              seq(0.2, 0.8, by = 0.05),
              seq(0.8, 1, by = 0.01))
    yell <- sqrt(xell * (1-xell))
    
    lab <- extract_label(node)
    
    # Fill node box color based on covariates type
    cov_type <- attr(meta[[sub("\\..*", "", lab[1])]], 'cov_type')
    fill[1] <- if(!is.null(cov_type)){
      switch(cov_type,
             fdata   = 'khaki1',
             diagram = 'wheat',
             graph   = 'darkseagreen1')
    } else {
      switch(class(meta[[sub("\\..*", "", lab[1])]]),
             logical = 'mistyrose',
             integer = 'lightsalmon',
             numeric = 'lightsalmon',
             factor  = 'aliceblue')
    }
    
    # Fill node id number box with a gray lighter than lightgray
    fill[2] <- 'gray96'
    
    grid.polygon(x = unit(c(xell, rev(xell)), "npc"),
                 y = unit(c(yell, -yell)+0.5, "npc"),
                 gp = gpar(fill = fill[1]))
    
    ## FIXME: something more general instead of pval ?
    grid.text(lab[1L], y = unit(1.5 + 0.5 * (lab[2L] != ""), "lines"))
    if(lab[2L] != "") grid.text(lab[2L], y = unit(1, "lines"))
    
    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
                           width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
                           height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2]))
      grid.text(nam[id_node(node)])
      popViewport()
    }
    upViewport(2)
  }
  
  return(rval)
}
class(node_inner) <- "grapcon_generator"


.nobs_party <- function(party, id = 1L) {
  dat <- data_party(party, id = id)
  if("(weights)" %in% names(dat)) sum(dat[["(weights)"]]) else NROW(dat)
}

edge_simple <- function(obj, digits = 3, abbreviate = FALSE, justmin = Inf,
                        just = c("alternate", "increasing", "decreasing",
                                 "equal"),
                        fill = "white")
{
  meta <- obj$data
  
  justfun <- function(i, split) {
    myjust <- if(mean(nchar(split)) > justmin) {
      match.arg(just, c("alternate", "increasing", "decreasing", "equal"))
    } else {
      "equal"
    }
    k <- length(split)
    rval <- switch(myjust,
                   "equal" = rep.int(0, k),
                   "alternate" = rep(c(0.5, -0.5), length.out = k),
                   "increasing" = seq(from = -k/2, to =  k/2, by = 1),
                   "decreasing" = seq(from =  k/2, to = -k/2, by = -1)
    )
    unit(0.5, "npc") + unit(rval[i], "lines")
  }
  
  ### panel function for simple edge labelling
  function(node, i) {
    #attempt to retrieve nobs for the kidnodes, instead of using regmatches etc:
    #print(nodeapply(obj, ids = id_node(node), FUN = function(n) nrow(n$data), by_node = FALSE))
    #would only take to be able to retrieve kidnodes' ids, instead of id_node(node)
    curr_split_type <- attr(split_node(node), 'curr_split_type')
    split_couple <- character_split(split_node(node), meta, digits = digits)$levels
    y <- justfun(i, split_couple)
    split <- split_couple[i]
    
    if (any(grep(">", split) > 0) || any(grep("<", split) > 0)) {
      # coeff or numeric/integer -> edge: parsed split
      
      #parse to turn '<=' into the symbol 'less than or equal to'
      #try() because parse(...) won't work for split = "< 10 Euro", for ex.
      tr <- suppressWarnings(try(parse(text = paste("phantom(0)", split)), silent = TRUE))
      #FIXME: parse adds unnecessary space to the left
      if(!inherits(tr, "try-error")) split <- tr
      
      grid.rect(y = y, gp = gpar(fill = fill, col = 0), width = unit(1.5, "strwidth", split))
      grid.text(split, y = y, just = "center")
      
    } else if (!is.null(curr_split_type) && curr_split_type == 'cluster'){
      #cluster -> edge: nobs in the kidnode
      
      #the number of obs in each kid node is calculated as the number of
      #commas appearing in split (which is a string where the levels are
      #separated by commas), plus one
      n_kid <- as.character(lengths(regmatches(split, gregexpr(",", split))) + 1)
      n_kid <- paste('n =', n_kid)
      
      grid.rect(y = y, gp = gpar(fill = fill, col = 0), width = unit(1, "strwidth", n_kid))
      grid.text(n_kid, y = y, just = "center")
      
    } else {
      #all the others -> edge: split
      
      grid.rect(y = y, gp = gpar(fill = fill, col = 0), width = unit(1, "strwidth", split))
      grid.text(split, y = y, just = "center")
      
    }
  }
}
class(edge_simple) <- "grapcon_generator"

.plot_node <- function(node, xlim, ylim, nx, ny,
                       terminal_panel, inner_panel, edge_panel,
                       tnex = 2, drop_terminal = TRUE, debug = FALSE) {
  
  ### the workhorse for plotting trees
  
  ### set up viewport for terminal node
  if (is.terminal(node)) {
    x <- xlim[1] + diff(xlim)/2
    y <- ylim[1] + 0.5
    
    tn_vp <- viewport(x = unit(x, "native"),
                      y = unit(y, "native") - unit(0.5, "lines"),
                      width = unit(1, "native"),
                      height = unit(tnex, "native") - unit(1, "lines"),
                      just = c("center", "top"),
                      name = paste("Node", id_node(node), sep = ""))
    pushViewport(tn_vp)
    if (debug)
      grid.rect(gp = gpar(lty = "dotted", col = 4))
    terminal_panel(node)
    upViewport()
    return(NULL)
  }
  
  ## convenience function for computing relative position of splitting node
  pos_frac <- function(node) {
    if(is.terminal(node)) 0.5 else {
      width_kids <- sapply(kids_node(node), width)
      nk <- length(width_kids)
      rval <- if(nk %% 2 == 0) sum(width_kids[1:(nk/2)]) else
        mean(cumsum(width_kids)[nk/2 + c(-0.5, 0.5)])
      rval/sum(width_kids)
    }
  }
  
  ## extract information
  split <- split_node(node)
  kids <- kids_node(node)
  width_kids <- sapply(kids, width)
  nk <- length(width_kids)
  
  ### position of inner node
  x0 <- xlim[1] + pos_frac(node) * diff(xlim)
  y0 <- max(ylim)
  
  ### relative positions of kids
  xfrac <- sapply(kids, pos_frac)
  x1lim <- xlim[1] + cumsum(c(0, width_kids))/sum(width_kids) * diff(xlim)
  x1 <- x1lim[1:nk] + xfrac * diff(x1lim)
  if (!drop_terminal) {
    y1 <- rep(y0 - 1, nk)
  } else {
    y1 <- ifelse(sapply(kids, is.terminal), tnex - 0.5, y0 - 1)
  }
  
  ### draw edges
  for(i in 1:nk) grid.lines(x = unit(c(x0, x1[i]), "native"), y = unit(c(y0, y1[i]), "native"))
  
  ### create viewport for inner node
  in_vp <- viewport(x = unit(x0, "native"),
                    y = unit(y0, "native"),
                    width = unit(1, "native"),
                    height = unit(1, "native") - unit(1, "lines"),
                    name = paste("Node", id_node(node), sep = ""))
  pushViewport(in_vp)
  if(debug) grid.rect(gp = gpar(lty = "dotted"))
  inner_panel(node)
  upViewport()
  
  ### position of labels
  y1max <- max(y1)
  ypos <- y0 - (y0 - y1max) * 0.5
  xpos <- x0 - (x0 - x1) * 0.5 * (y0 - y1max)/(y0 - y1)
  
  ### setup labels
  for(i in 1:nk) {
    sp_vp <- viewport(x = unit(xpos[i], "native"),
                      y = unit(ypos, "native"),
                      width = unit(diff(x1lim)[i], "native"),
                      height = unit(1, "lines"),
                      name =  paste("edge", id_node(node), "-", i, sep = ""))
    pushViewport(sp_vp)
    if(debug) grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_panel(node, i)
    upViewport()
  }
  
  ## call workhorse for kids
  for(i in 1:nk) .plot_node(kids[[i]],
                            c(x1lim[i], x1lim[i+1]), c(y1[i], 1), nx, ny,
                            terminal_panel, inner_panel, edge_panel,
                            tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}

node_barplot <- function(obj,
                         col = "black",
                         fill = NULL,
                         bg = "white",
                         beside = NULL,
                         ymax = NULL,
                         ylines = NULL,
                         widths = 1,
                         gap = NULL,
                         reverse = NULL,
                         rot = 0,
                         just = c("center", "top"),
                         id = TRUE,
                         mainlab = NULL,
                         text = c("none", "horizontal", "vertical"),
                         gp = gpar())
{
  ## extract response
  y <- obj$fitted[["(response)"]]
  stopifnot(is.factor(y) || isTRUE(all.equal(round(y), y)) || is.data.frame(y))
  
  ## FIXME: This could be avoided by
  ##   predict_party(obj, nodeids(obj, terminal = TRUE), type = "prob")
  ## but only for terminal nodes                  ^^^^
  probs_and_n <- function(x) {
    y1 <- x$fitted[["(response)"]]
    if(!is.factor(y1)) {
      if(is.data.frame(y1)) {
        y1 <- t(as.matrix(y1))
      } else {
        y1 <- factor(y1, levels = min(y):max(y))
      }
    }
    w <- x$fitted[["(weights)"]]
    if(is.null(w)) w <- rep.int(1L, length(y1))
    sumw <- if(is.factor(y1)) tapply(w, y1, sum) else drop(y1 %*% w)
    sumw[is.na(sumw)] <- 0
    prob <- c(sumw/sum(w), sum(w))
    names(prob) <- c(if(is.factor(y1)) levels(y1) else rownames(y1), "nobs")
    prob
  }
  probs <- do.call("rbind", nodeapply(obj, nodeids(obj), probs_and_n, by_node = FALSE))
  nobs <- probs[, "nobs"]
  probs <- probs[, -ncol(probs), drop = FALSE]
  if(is.factor(y)) {
    ylevels <- levels(y)
    if(is.null(beside)) beside <- if(length(ylevels) < 3L) FALSE else TRUE
    if(is.null(ymax)) ymax <- if(beside) 1.1 else 1
    if(is.null(gap)) gap <- if(beside) 0.1 else 0
  } else {
    if(is.null(beside)) beside <- TRUE
    if(is.null(ymax)) ymax <- if(beside) max(probs) * 1.1 else max(probs)
    ylevels <- colnames(probs)
    if(length(ylevels) < 2) ylevels <- ""
    if(is.null(gap)) gap <- if(beside) 0.1 else 0
  }
  if(is.null(reverse)) reverse <- !beside
  if(is.null(fill)) fill <- gray.colors(length(ylevels))
  if(is.null(ylines)) ylines <- if(beside) c(3, 2) else c(1.5, 2.5)
  
  ## text labels?
  if(isTRUE(text)) text <- "horizontal"
  if(!is.character(text)) text <- "none"
  text <- match.arg(text, c("none", "horizontal", "vertical"))
  
  ### panel function for barplots in nodes
  rval <- function(node) {
    
    ## id
    nid <- id_node(node)
    
    ## parameter setup
    pred <- probs[nid,]
    if(reverse) {
      pred <- rev(pred)
      ylevels <- rev(ylevels)
    }
    np <- length(pred)
    nc <- if(beside) np else 1
    
    fill <- rep(fill, length.out = np)
    widths <- rep(widths, length.out = nc)
    col <- rep(col, length.out = nc)
    ylines <- rep(ylines, length.out = 2)
    
    gap <- gap * sum(widths)
    yscale <- c(0, ymax)
    xscale <- c(0, sum(widths) + (nc+1)*gap)
    
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste0("node_barplot", nid),
                       gp = gp)
    
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    
    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(names(obj)[nid], nobs[nid])
    }
    grid.text(mainlab)
    popViewport()
    
    plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                     xscale=xscale, yscale=yscale,
                     name = paste0("node_barplot", node$nodeID, "plot"),
                     clip = FALSE)
    
    pushViewport(plot)
    
    if(beside) {
      xcenter <- cumsum(widths+gap) - widths/2
      if(length(xcenter) > 1) grid.xaxis(at = xcenter, label = FALSE)
      grid.text(ylevels, x = xcenter, y = unit(-1, "lines"),
                just = just, rot = rot,
                default.units = "native", check.overlap = TRUE)
      grid.yaxis()
      grid.rect(gp = gpar(fill = "transparent"))
      grid.clip()
      for (i in 1:np) {
        grid.rect(x = xcenter[i], y = 0, height = pred[i],
                  width = widths[i],
                  just = c("center", "bottom"), default.units = "native",
                  gp = gpar(col = col[i], fill = fill[i]))
        if(text != "none") {
          grid.text(x = xcenter[i], y = pred[i] + 0.025,
                    label = paste(format(round(100 * pred[i], 1), nsmall = 1), "%", sep = ""),
                    just = if(text == "horizontal") c("center", "bottom") else c("left", "center"),
                    rot = if(text == "horizontal") 0 else 90,
                    default.units = "native")
        }
      }
    } else {
      ycenter <- cumsum(pred) - pred
      
      if(np > 1) {
        grid.text(ylevels[1], x = unit(-1, "lines"), y = 0,
                  just = c("left", "center"), rot = 90,
                  default.units = "native", check.overlap = TRUE)
        grid.text(ylevels[np], x = unit(-1, "lines"), y = ymax,
                  just = c("right", "center"), rot = 90,
                  default.units = "native", check.overlap = TRUE)
      }
      if(np > 2) {
        grid.text(ylevels[-c(1,np)], x = unit(-1, "lines"), y = ycenter[-c(1,np)],
                  just = "center", rot = 90,
                  default.units = "native", check.overlap = TRUE)
      }
      grid.yaxis(main = FALSE)
      
      grid.clip()
      grid.rect(gp = gpar(fill = "transparent"))
      for (i in 1:np) {
        grid.rect(x = xscale[2]/2, y = ycenter[i], height = min(pred[i], ymax - ycenter[i]),
                  width = widths[1],
                  just = c("center", "bottom"), default.units = "native",
                  gp = gpar(col = col[i], fill = fill[i]))
      }
    }
    grid.rect(gp = gpar(fill = "transparent"))
    
    
    upViewport(2)
  }
  
  return(rval)
}
class(node_barplot) <- "grapcon_generator"

node_boxplot <- function(obj,
                         col = "black",
                         fill = "lightgray",
                         bg = "white",
                         width = 0.5,
                         yscale = NULL,
                         ylines = 3,
                         cex = 0.5,
                         id = TRUE,
                         mainlab = NULL,
                         gp = gpar())
{
  y <- obj$fitted[["(response)"]]
  stopifnot(is.numeric(y))
  
  if (is.null(yscale))
    yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
  
  ### panel function for boxplots in nodes
  rval <- function(node) {
    
    ## extract data
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    yn <- dat[["(response)"]]
    wn <- dat[["(weights)"]]
    if(is.null(wn)) wn <- rep(1, length(yn))
    
    ## parameter setup
    x <- boxplot(rep.int(yn, wn), plot = FALSE)
    
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1),
                                                          c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_boxplot", nid, sep = ""),
                       gp = gp)
    
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    
    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(names(obj)[nid], sum(wn))
    }
    grid.text(mainlab)
    popViewport()
    
    plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                     xscale = c(0, 1), yscale = yscale,
                     name = paste0("node_boxplot", nid, "plot"),
                     clip = FALSE)
    
    pushViewport(plot)
    
    grid.yaxis()
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()
    
    xl <- 0.5 - width/4
    xr <- 0.5 + width/4
    
    ## box & whiskers
    grid.lines(unit(c(xl, xr), "npc"),
               unit(x$stats[1], "native"), gp = gpar(col = col))
    grid.lines(unit(0.5, "npc"),
               unit(x$stats[1:2], "native"), gp = gpar(col = col, lty = 2))
    grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"),
              width = unit(width, "npc"), height = unit(diff(x$stats[c(2, 4)]), "native"),
              just = c("center", "bottom"),
              gp = gpar(col = col, fill = fill))
    grid.lines(unit(c(0.5 - width/2, 0.5+width/2), "npc"),
               unit(x$stats[3], "native"), gp = gpar(col = col, lwd = 2))
    grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"),
               gp = gpar(col = col, lty = 2))
    grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"),
               gp = gpar(col = col))
    
    ## outlier
    n <- length(x$out)
    if (n > 0) {
      index <- 1:n ## which(x$out > yscale[1] & x$out < yscale[2])
      if (length(index) > 0)
        grid.points(unit(rep.int(0.5, length(index)), "npc"),
                    unit(x$out[index], "native"),
                    size = unit(cex, "char"), gp = gpar(col = col))
    }
    
    upViewport(2)
  }
  
  return(rval)
}
class(node_boxplot) <- "grapcon_generator"

node_surv <- function(obj, col = "black", bg = "white", yscale = c(0, 1), ylines = 2,
                      id = TRUE, mainlab = NULL, gp = gpar(), ...)
{
  ## extract response
  y <- obj$fitted[["(response)"]]
  stopifnot(inherits(y, "Surv"))
  
  ## helper functions
  mysurvfit <- function(y, weights, ...)
    survfit(y ~ 1, weights = weights)
  ### structure(
  ###   survival:::survfitKM(x = gl(1, NROW(y)), y = y, casewt = weights, ...),
  ### class = "survfit")
  
  dostep <- function(x, y) {
    ### create a step function based on x, y coordinates
    ### modified from `survival:print.survfit'
    if (is.na(x[1] + y[1])) {
      x <- x[-1]
      y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {
      # replace verbose horizonal sequences like
      # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
      # with (1, .2), (3, .1).  They are slow, and can smear the looks
      # of the line type.
      dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
      n2 <- sum(dupy)
      
      #create a step function
      xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
      yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
      RET <- list(x = xrep, y = yrep)
    } else {
      if (n == 1) {
        RET <- list(x = x, y = y)
      } else {
        RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
      }
    }
    return(RET)
  }
  
  ### panel function for Kaplan-Meier curves in nodes
  rval <- function(node) {
    
    ## extract data
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    yn <- dat[["(response)"]]
    wn <- dat[["(weights)"]]
    if(is.null(wn)) wn <- rep(1, NROW(yn))
    
    ## get Kaplan-Meier curver in node
    km <- mysurvfit(yn, weights = wn, ...)
    a <- dostep(km$time, km$surv)
    
    ## set up plot
    yscale <- yscale
    xscale <- c(0, max(y[,1]))
    
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1),
                                                          c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_surv", nid, sep = ""), gp = gp)
    
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    
    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, sum(wn))
    }
    grid.text(mainlab)
    popViewport()
    
    plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                     xscale=xscale, yscale = yscale,
                     name = paste0("node_surv", nid, "plot"),
                     clip = FALSE)
    
    pushViewport(plot)
    grid.xaxis()
    grid.yaxis()
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()
    grid.lines(unit(a$x, "native"), unit(a$y, "native"), gp = gpar(col = col))
    upViewport(2)
  }
  
  return(rval)
}
class(node_surv) <- "grapcon_generator"

node_ecdf <- function(obj, col = "black", bg = "white", ylines = 2,
                      id = TRUE, mainlab = NULL, gp = gpar(), ...)
{
  ## extract response
  y <- obj$fitted[["(response)"]]
  stopifnot(inherits(y, "numeric") || inherits(y, "integer"))
  
  dostep <- function(f) {
    x <- knots(f)
    y <- f(x)
    ### create a step function based on x, y coordinates
    ### modified from `survival:print.survfit'
    if (is.na(x[1] + y[1])) {
      x <- x[-1]
      y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {
      # replace verbose horizonal sequences like
      # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
      # with (1, .2), (3, .1).  They are slow, and can smear the looks
      # of the line type.
      dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
      n2 <- sum(dupy)
      
      #create a step function
      xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
      yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
      RET <- list(x = xrep, y = yrep)
    } else {
      if (n == 1) {
        RET <- list(x = x, y = y)
      } else {
        RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
      }
    }
    return(RET)
  }
  
  ### panel function for ecdf in nodes
  rval <- function(node) {
    
    ## extract data
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    yn <- dat[["(response)"]]
    wn <- dat[["(weights)"]]
    if(is.null(wn)) wn <- rep(1, NROW(yn))
    
    ## get ecdf in node
    f <- .pred_ecdf(yn, wn)
    a <- dostep(f)
    
    ## set up plot
    yscale <- c(0, 1)
    xscale <- range(y)
    a$x <- c(xscale[1], a$x[1], a$x, xscale[2])
    a$x <- a$x - min(a$x)
    a$x <- a$x / max(a$x)
    a$y <- c(0, 0, a$y, 1)
    
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1),
                                                          c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_ecdf", nid, sep = ""), gp = gp)
    
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    
    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, sum(wn))
    }
    grid.text(mainlab)
    popViewport()
    
    plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                     xscale=xscale, yscale=yscale,
                     name = paste0("node_surv", nid, "plot"),
                     clip = FALSE)
    
    pushViewport(plot)
    grid.xaxis()
    grid.yaxis()
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()
    grid.lines(a$x, a$y, gp = gpar(col = col))
    upViewport(2)
  }
  
  return(rval)
}
class(node_ecdf) <- "grapcon_generator"



node_mvar <- function(obj, which = NULL, id = TRUE, pop = TRUE, ylines = NULL, mainlab = NULL, varlab = TRUE, bg = "white", ...)
{
  ## obtain dependent variables
  y <- obj$fitted[["(response)"]]
  
  ## fitted node ids
  fitted <- obj$fitted[["(fitted)"]]
  
  ## number of panels needed
  if(is.null(which)) which <- 1L:NCOL(y)
  k <- length(which)
  
  rval <- function(node) {
    
    tid <- id_node(node)
    nobs <- .nobs_party(obj, id = tid)
    
    ## set up top viewport
    top_vp <- viewport(layout = grid.layout(nrow = k, ncol = 2,
                                            widths = unit(c(ylines, 1), c("lines", "null")), heights = unit(k, "null")),
                       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_mvar", tid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    
    ## main title
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(tid, nobs)
    }
    
    for(i in 1L:k) {
      tmp <- obj
      tmp$fitted[["(response)"]] <- y[,which[i]]
      if(varlab) {
        nm <- names(y)[which[i]]
        if(i == 1L) nm <- paste(mainlab, nm, sep = ": ")
      } else {
        nm <- if(i == 1L) mainlab else ""
      }
      pfun <- switch(sapply(y, class)[which[i]],
                     "Surv" = node_surv(tmp, id = id, mainlab = nm, ...),
                     "factor" = node_barplot(tmp, id = id, mainlab = nm,  ...),
                     "ordered" = node_barplot(tmp, id = id, mainlab = nm, ...),
                     node_boxplot(tmp, id = id, mainlab = nm, ...))
      ## select panel
      plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = i)
      pushViewport(plot_vpi)
      
      ## call panel function
      pfun(node)
      
      if(pop) popViewport() else upViewport()
    }
    if(pop) popViewport() else upViewport()
  }
  
  return(rval)
}
class(node_mvar) <- "grapcon_generator"


# split -------------------------------------------------------------------

partysplit <- function(varid, breaks = NULL, index = NULL, right = TRUE,
                       prob = NULL, info = NULL, centroids = NULL, basid = NULL) {
  
  ### informal class for splits
  split <- vector(mode = "list", length = 8)
  names(split) <- c("varid", "breaks", "index", "right", "prob", "info", "centroids", "basid")
  
  ### split is an id referring to a variable
  if (!is.integer(varid))
    stop(sQuote("varid"), " ", "is not integer")
  split$varid <- varid
  
  if (is.null(breaks) && is.null(index))
    stop("either", " ", sQuote("breaks"), " ", "or", " ",
         sQuote("index"), " ", "must be given")
  
  ### vec
  if (!is.null(breaks)) {
    if (is.numeric(breaks) && (length(breaks) >= 1)) {
      ### FIXME: I think we need to make sure breaks are double in C
      split$breaks <- as.double(breaks)
    } else {
      stop(sQuote("break"), " ",
           "should be a numeric vector containing at least one element")
    }
  }

  if (!is.null(index)) {
    if (is.integer(index)) {
      if (!(length(index) >= 2))
        stop(sQuote("index"), " ", "has less than two elements")
      if (!(min(index, na.rm = TRUE) == 1))
        stop("minimum of", " ", sQuote("index"), " ", "is not equal to 1")
      if (!all.equal(diff(sort(unique(index))), rep(1, max(index, na.rm = TRUE) - 1)))
        stop(sQuote("index"), " ", "is not a contiguous sequence")
      split$index <- index
    } else {
      stop(sQuote("index"), " ", "is not a class", " ", sQuote("integer"))
    }
    if (!is.null(breaks)) {
      if (length(breaks) != (length(index) - 1))
        stop("length of", " ", sQuote("breaks"), " ",
             "does not match length of", " ", sQuote("index"))
    }
  }
  
  if (is.logical(right) & !is.na(right))
    split$right <- right
  else
    stop(sQuote("right"), " ", "is not a logical")
  
  if (!is.null(prob)) {
    if (!is.double(prob) ||
        (any(prob < 0) | any(prob > 1) | !isTRUE(all.equal(sum(prob), 1))))
      stop(sQuote("prob"), " ", "is not a vector of probabilities")
    if (!is.null(index)) {
      if (!(max(index, na.rm = TRUE) == length(prob)))
        stop("incorrect", " ", sQuote("index"))
    }
    if (!is.null(breaks) && is.null(index)) {
      if (!(length(breaks) == (length(prob) - 1)))
        stop("incorrect", " ", sQuote("breaks"))
    }
    split$prob <- prob
  }
  
  if (!is.null(info))
    split$info <- info
  
  if (!is.null(centroids))
    split$centroids <- centroids
  
  if (!is.null(basid))
    if (!is.integer(varid)){
      stop(sQuote("varid"), " ", "is not integer")
    } else {
      split$basid <- basid
    }
  
  class(split) <- "partysplit"
  
  return(split)
}

centroids_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$centroids
}

basid_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$basid
}


kidids_split <- function(split, data, vmatch = 1:length(data), obs = NULL) {
  
  varid <- varid_split(split)
  basid <- basid_split(split)
  class(data) <- "list" ### speed up
  
  if(!is.null(basid)){ #coeff
    x <- data[[vmatch[varid]]][[as.character(basid)]]
  } else {
    x <- data[[vmatch[varid]]]
  }
  if (!is.null(obs)) x <- x[obs]
  
  if (is.null(breaks_split(split))) {
    if (storage.mode(x) != "integer")
      stop("variable", " ", vmatch[varid], " ", "is not integer")
  } else {
    ### labels = FALSE returns integers and is faster
    ### <FIXME> use findInterval instead of cut?
    #        x <- cut.default(as.numeric(x), labels = FALSE,
    #                 breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
    #                 right = right_split(split))
    x <- .bincode(as.numeric(x), # labels = FALSE,
                  breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
                  right = right_split(split))
    ### </FIXME>
  }
  index <- index_split(split)
  ### empty factor levels correspond to NA and return NA here
  ### and thus the corresponding observations will be treated
  ### as missing values (surrogate or random splits):
  if (!is.null(index))
    x <- index[x]
  return(x)
}

kidids_split_predict <- function(split, data, vmatch = 1:length(data), obs = NULL) {
  
  varid <- varid_split(split)
  basid <- basid_split(split)
  class(data) <- "list" ### speed up
  
  if(!is.null(basid)){ #coeff
    x <- data[[vmatch[varid]]][[as.character(basid)]]
  } else {
    x <- data[[vmatch[varid]]]
  }
  if (!is.null(obs)) x <- x[obs]
  
  if (is.null(breaks_split(split))) {
    
    if(!is.null(centroids_split(split))){
      
      centroids <- centroids_split(split)
      centroids_dist <- sapply(centroids, dist_comp_cl, x = x, lp = 2)
      x <- if(!is.null(dim(centroids_dist))) apply(centroids_dist, 1, which.min)
      else which.min(centroids_dist)
      
    }
    
    if (storage.mode(x) != "integer")
      stop("variable", " ", vmatch[varid], " ", "is not integer")
  } else {
    ### labels = FALSE returns integers and is faster
    ### <FIXME> use findInterval instead of cut?
    #        x <- cut.default(as.numeric(x), labels = FALSE,
    #                 breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
    #                 right = right_split(split))
    x <- .bincode(as.numeric(x), # labels = FALSE,
                  breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
                  right = right_split(split))
    ### </FIXME>
  }
  index <- index_split(split)
  ### empty factor levels correspond to NA and return NA here
  ### and thus the corresponding observations will be treated
  ### as missing values (surrogate or random splits):
  if (!is.null(index) & is.null(centroids_split(split)))
    x <- index[x]
  return(x)
}

character_split <- function(split, data = NULL, digits = getOption("digits") - 2) {
  
  varid <- varid_split(split)
  basid <- basid_split(split)

  if (!is.null(data)) {
    ## names and labels
    lev <- lapply(data, levels)[[varid]]
    mlab <- names(data)[varid]

    ## determine split type
    type <- sapply(data, function(x) class(x)[1])[varid_split(split)]
    type[!(type %in% c("factor", "ordered"))] <- "numeric"
  } else {
    ## (bad) default names and labels
    lev <- NULL
    mlab <- paste("V", varid, sep = "")
    type <- "numeric"
  }
  
  ## process defaults for breaks and index
  breaks <- breaks_split(split)
  index <- index_split(split)
  right <- right_split(split)
  
  if (is.null(breaks)) breaks <- 1:(length(index) - 1)
  if (is.null(index)) index <- 1:(length(breaks) + 1)
  
  ## check whether ordered are really ordered
  if (type == "ordered") {
    if (length(breaks) > 1)
      type <- "factor"
  }
  ### <FIXME> format ordered multiway splits? </FIXME>
  
  switch(type, 
         "factor" = {
           nindex <- index[cut(seq_along(lev), c(-Inf, breaks, Inf), right = right)]
           dlab <- as.vector(tapply(lev, nindex, paste, collapse = ", "))
         },
         "ordered" = {
           if (length(breaks) == 1) {
             if (right)
               dlab <- paste(c("<=", ">"), lev[breaks], sep = " ")
             else
               dlab <- paste(c("<", ">="), lev[breaks], sep = " ")
           } else {
             stop("") ### see above
           }
           dlab <- dlab[index]
         },
         "numeric" = {
           breaks <- round(breaks, digits)
           if (length(breaks) == 1) {
             if (right)
               dlab <- paste(c("<=", ">"), breaks, sep = " ")
             else
               dlab <- paste(c("<", ">="), breaks, sep = " ")
           } else {
             dlab <- levels(cut(0, breaks = c(-Inf, breaks, Inf), 
                                right = right))
           }
           dlab <- as.vector(tapply(dlab, index, paste, collapse = " | "))
         }
  )

  return(list(name = mlab, basid = basid, levels = dlab))
}


# print -------------------------------------------------------------------

print.partynode <- function(x, data = NULL, names = NULL,
                            inner_panel = function(node) "", terminal_panel = function(node) " *",
                            prefix = "", first = TRUE, digits = getOption("digits") - 2, ...)
{
  ids <- nodeids(x)
  
  if(first) {
    if(is.null(names)) names <- as.character(ids)
    cat(paste(prefix, "[", names[which(ids == id_node(x))], "] root", sep = ""))
    
    if(is.terminal(x)) {
      char <- terminal_panel(x)
      if(length(char) > 1L) {
        cat(paste(char[1L], "\n",
                  paste(prefix, "    ", char[-1L], sep = "", collapse = "\n"),
                  sep = ""), "\n")
      } else {
        cat(char, "\n")
      }
    } else {
      cat("\n")
    }
  }
  
  if (length(x) > 0) {
    ## add indentation
    nextprefix <- paste(prefix, "|   ", sep = "")
    
    ## split labels
    slabs <- character_split(split_node(x), data = data, digits = digits, ...)
    slabs <- ifelse(substr(slabs$levels, 1, 1) %in% c("<", ">"),
                    paste(if(!is.null(slabs$basid))
                      paste0(slabs$name, '.', slabs$basid) else slabs$name,
                      slabs$levels), 
                    paste(if(!is.null(slabs$basid))
                      paste0(slabs$name, '.', slabs$basid) else slabs$name,
                      "in", slabs$levels))
    
    ## kid labels
    knodes <- kids_node(x)    
    knam <- sapply(knodes, function(z) names[which(ids == id_node(z))])
    klabs <- sapply(knodes, function(z)
      if(is.terminal(z)) {
        char <- terminal_panel(z)
        if(length(char) > 1L) {
          paste(char[1L], "\n",
                paste(nextprefix, "    ", char[-1L], sep = "", collapse = "\n"),
                sep = "")
        } else {
          char
        }
      } else {
        paste(inner_panel(z), collapse = "\n")
      })
    
    ## merge, cat, and call recursively
    labs <- paste("|   ", prefix, "[", knam, "] ", slabs, klabs, "\n", sep = "")          
    for (i in 1:length(x)) {
      cat(labs[i])
      print.partynode(x[i], data = data, names = names[match(nodeids(x[i]), ids)],
                      inner_panel = inner_panel, terminal_panel = terminal_panel,
                      prefix = nextprefix,  first = FALSE, digits = digits, ...)
    }
  }
}

print.party <- function(x,
                        terminal_panel = function(node) formatinfo_node(node, default = "*", prefix = ": "), tp_args = list(),
                        inner_panel = function(node) "", ip_args = list(),
                        header_panel = function(party) "",
                        footer_panel = function(party) "",
                        digits = getOption("digits") - 2, ...)
{
  ## header
  cat(paste(header_panel(x), collapse = "\n"))

  ## nodes
  if(inherits(terminal_panel, "grapcon_generator"))
    terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
  if(inherits(inner_panel, "grapcon_generator"))
    inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
  print(node_party(x), x$data, names = names(x),
        terminal_panel = terminal_panel, inner_panel = inner_panel,
        digits = digits, ...)
  
  ## footer
  cat(paste(footer_panel(x), collapse = "\n"))
}

print.simpleparty <- function(x, digits = getOption("digits") - 4,
                              header = NULL, footer = TRUE, ...)
{
  ## header panel
  if(is.null(header)) header <- !is.null(terms(x))
  header_panel <- if(header) function(party) {
    c("", "Model formula:", deparse(formula(terms(party))), "", "Fitted party:", "")
  } else function(party) ""
  
  ## footer panel
  footer_panel <- if(footer) function(party) {
    n <- width(party)
    n <- format(c(length(party) - n, n))
    
    c("", paste("Number of inner nodes:   ", n[1]),
      paste("Number of terminal nodes:", n[2]), "")
  } else function (party) ""
  
  ## terminal panel
  terminal_panel <- function(node) formatinfo_node(node,
                                                   FUN = .make_formatinfo_simpleparty(x, digits = digits),
                                                   default = "*", prefix = ": ")
  
  print.party(x, terminal_panel = terminal_panel,
              header_panel = header_panel, footer_panel = footer_panel, ...)
}


# other -------------------------------------------------------------------

.make_formatinfo_simpleparty <- function(x, digits = getOption("digits") - 4, sep = "")
{
  ## digit processing
  digits <- max(c(0, digits))
  digits2 <- max(c(0, digits - 2))
  
  ## type of predictions
  y <- node_party(x)$info$prediction
  yclass <- class(y)[1]
  if(yclass == "ordered") yclass <- "factor"
  if(!(yclass %in% c("survfit", "factor"))) yclass <- "numeric"
  
  ## type of weights
  n <- node_party(x)$info$n
  if(is.null(names(n))) {
    wdigits <- 0
    wsym <- "n"
  } else {
    if(names(n) == "w") {
      wdigits <- max(c(0, digits - 2))
      wsym <- "w"
    } else {
      wdigits <- 0
      wsym <- "n"
    }
  }
  
  ## compute terminal node labels
  FUN <- function(info) {
    yhat <- info$prediction
    if (yclass == "survfit") {
      yhat <- .median_survival_time(yhat)
      yclass <- "numeric"
    }
    if(yclass == "numeric") yhat <- format(round(yhat, digits = digits), nsmall = digits)
    w <- info$n
    yerr <- if(is.null(info$error)) "" else paste(", err = ",
                                                  format(round(info$error, digits = digits2), nsmall = digits2),
                                                  names(info$error), sep = "")
    rval <- paste(yhat, sep,
                  " (", wsym, " = ", format(round(w, digits = wdigits), nsmall = wdigits),
                  yerr, ")", sep = "")
    unlist(strsplit(rval, "\n"))
  }
  return(FUN)
}

.resample <- function(x, ...) x[sample.int(length(x), ...)]

.median_survival_time <- function(x) {
  minmin <- function(y, xx) {
    if (any(!is.na(y) & y==.5)) {
      if (any(!is.na(y) & y <.5))
        .5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
      else
        .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
    } else   min(xx[!is.na(y) & y<=.5])
  }
  med <- suppressWarnings(minmin(x$surv, x$time))
  return(med)
}

### get the recursive index
### obj is of class "partynode"
.get_path <- function(obj, i) {
  
  idx <- c()
  recFun <- function(node, i) {
    if (id_node(node) == i) return(NULL)
    idx <<- c(idx, which(names(unclass(node)) == "kids"))
    kid <- sapply(kids_node(node), id_node)
    nextid <- max(which(kid <= i))
    idx <<- c(idx, nextid)
    return(recFun(node[[nextid]], i))
  }
  out <- recFun(obj, i)
  return(idx)
}
