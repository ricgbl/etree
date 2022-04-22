findsplit <- function(response,
                      covariates,
                      alpha,
                      R,
                      split_type,
                      coeff_split_type,
                      p_adjust_method,
                      random_covs) {
  
  # Subset of covariates to consider for splitting (possibly all of them)
  cov_subset <- as.integer(seq_along(covariates$cov))
  if (!is.null(random_covs)) {
    cov_subset <- sample(cov_subset, random_covs)
  }
  
  # Independence test between the response and each covariate
  resp_dist <- response$response_dist
  stat_pval <- sapply(covariates$dist[cov_subset],
                      function(cov_dist) {
                        indep_test(x_dist = cov_dist, y_dist = resp_dist, R = R)
                      }
  )
  if (all(is.na(stat_pval['Pvalue', ]))) {
    return(NULL)
  }
  
  # Multiple testing correction
  adj_p <- p.adjust(stat_pval['Pvalue', ], method = p_adjust_method)
  
  # Stop criterion
  if (min(adj_p, na.rm = TRUE) > alpha) {
    return(NULL)
  }
  
  # Variable selection (based on original pvalues)
  xselect <- .select_variable(statistic_pvalue = stat_pval)
  #.select_variable is just a copy of .select_component
  
  # Selected covariates
  xselect <- cov_subset[xselect] #useful only if !is.null(random_covs)
  x <- covariates$cov[[xselect]]
  newx <- covariates$newcov[[xselect]]
  if (split_type == 'cluster') {
    xdist <- covariates$dist[[xselect]]
  }
  
  # Control for newx not being void
  if ((split_type == 'coeff') &&
      (class(x) == 'list') && all(sapply(x, class) == 'igraph') &&
      (dim(newx)[2] == 0)) {
    
    # Throw warning: if x is selected (through distance) and newx is void because
    # it had all-equal columns (after coeff exp), it means that either distance
    # or coeff exp are not appropriate, or at least that they are uncompatible
    warning('The selected covariate has non-informative coefficient expansion.
    Please reconsider the choice for at least one between distance and
    coefficient expansion. The selected covariate is ignored for the time
            being.')
    
    # Ignore the selected covariate and re-run findsplit
    # in order to keep the original indices (i.e. not to create confusion), the
    # selected covariate is not dropped, but it is instead replaced with a
    # vector of 0s, so that it is selected no more
    covariates$dist[[xselect]][] <- 0L
    split <- findsplit(response = response,
                       covariates = covariates,
                       alpha = alpha,
                       R = R,
                       split_type = split_type,
                       coeff_split_type = coeff_split_type,
                       p_adjust_method = p_adjust_method,
                       random_covs = random_covs)
    return(split)
  }
  
  # Split point search
  split_objs <- split_opt(resp = response,
                          cov = x,
                          new_cov = newx,
                          cov_dist = xdist,
                          R = R,
                          split_type = split_type,
                          coeff_split_type = coeff_split_type)
  
  # If split_objs is VOID, ignore the selected covariate and re-run findsplit
  if (isTRUE(split_objs$void)) {
    covariates$dist[[xselect]][] <- 0L
    split <- findsplit(response = response,
                       covariates = covariates,
                       alpha = alpha,
                       R = R,
                       split_type = split_type,
                       coeff_split_type = coeff_split_type,
                       p_adjust_method = p_adjust_method,
                       random_covs = random_covs)
    return(split)
  }
  
  # Separately save split_objs outputs
  splitpoint <- split_objs$splitpoint
  splitindex <- split_objs$splitindex
  bselect <- split_objs$bselect
  centroids <- split_objs$centroids

  # Return the split point
  switch(class(x),
         
         integer = {
           
           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitpoint,
                                  info = list(pvalue = stat_pval),
                                  right = TRUE))
           
         },
         
         numeric = {
           
           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitpoint,
                                  info = list(pvalue = stat_pval),
                                  right = TRUE))
           
         },
         
         factor = {
           
           return(sp = partysplit(varid = as.integer(xselect),
                                  index = splitindex,
                                  info = list(pvalue = stat_pval)))
           
         },
         
         fdata = {
           
           if (split_type == 'coeff') {
             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitpoint,
                                    right = TRUE,
                                    info = list(pvalue = stat_pval)))
             
           } else if (split_type == 'cluster') {
             
             sp <- partysplit(varid = as.integer(xselect),
                              centroids = centroids,
                              index = as.integer(splitindex),
                              info = list(pvalue = stat_pval))
             attr(sp, 'curr_split_type') <- 'cluster'   #used in edge.simple
             return(sp)
             
           }
         },
         
         list = if (all(sapply(x, function(x) attributes(x)$names)
                        == 'diagram')) {
           
           #only cluster
           sp <- partysplit(varid = as.integer(xselect),
                            centroids = centroids,
                            index = as.integer(splitindex),
                            info = list(pvalue = stat_pval))
           attr(sp, 'curr_split_type') <- 'cluster'   #used in edge.simple
           return(sp)
           
         } else if (all(sapply(x, class) == 'igraph')) {
           
           if (split_type == 'coeff') {
             
             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitpoint,
                                    right = TRUE,
                                    info = list(pvalue = stat_pval)))
             
           } else if (split_type == 'cluster') {
             
             sp <- partysplit(varid = as.integer(xselect),
                              centroids = centroids,
                              index = as.integer(splitindex),
                              info = list(pvalue = stat_pval))
             attr(sp, 'curr_split_type') <- 'cluster'   #used in edge.simple
             return(sp)
             
           }
         })
}


split_opt <- function(resp,
                      cov,
                      new_cov,
                      cov_dist,
                      R = 1000,
                      split_type = 'coeff',
                      coeff_split_type = 'test') {
  
  # Retrieve response_dist and response
  resp_dist <- resp$response_dist
  resp <- resp$response
  
  switch(class(cov),
         
         integer    = {
           
           s  <- sort(cov)
           comb <- sapply(s[-length(s)], function(j) cov <= j)
           
           if (coeff_split_type == 'traditional') {
             
             splitpoint <- .select_traditional(values = s,
                                               comb_logical = comb,
                                               resp = resp)
             
           } else if (coeff_split_type == 'test') {
             
             stat_pval <- apply(comb, 2,
                                function(q) indep_test(x = q,
                                                       y_dist = resp_dist,
                                                       R = R))
             splitpoint <- .select_splitpoint(values = s,
                                              statistic_pvalue = stat_pval)
             
           }
         },
         
         numeric    = {
           
           s  <- sort(cov)
           comb <- sapply(s[-length(s)], function(j) cov <= j)
           
           if (coeff_split_type == 'traditional') {
             
             splitpoint <- .select_traditional(values = s,
                                               comb_logical = comb,
                                               resp = resp)
             
           } else if (coeff_split_type == 'test') {
             
             stat_pval <- apply(comb, 2,
                                function(q) indep_test(x = q,
                                                       y_dist = resp_dist,
                                                       R = R))
             splitpoint <- .select_splitpoint(values = s,
                                              statistic_pvalue = stat_pval)
             
           }
         },
         
         factor     = {
           
           # Drop unused levels
           lev <- levels(cov[drop = TRUE])
           n_lev <- length(lev)
           
           if (n_lev == 2) {
             
             splitpoint <- lev[1]  #split point simply given by the first level
             
           } else {
             
             # Combination of all the levels
             lev_cmb <- do.call("c",
                                lapply(seq_len(n_lev / 2), function(ntaken) {
                                  #1 for n_lev=2,3; 1,2 for n_lev=4,5; 1,2,3 for
                                  #n_lev=6,7;... since [(1:x) := (1:floor(x))]
                                  combs_ntaken <- combn(x = lev,
                                                        m = ntaken,
                                                        simplify = FALSE)
                                  #fine, when n_lev is odd; for even n_lev, skip
                                  #half the combs for the last ntaken
                                  #e.g. n_lev=6 => last_ntaken=3, but '1,5,6' is
                                  #equivalent to '2,3,4'
                                  if (ntaken == (n_lev / 2)) {
                                    combs_ntaken <-
                                      combs_ntaken[seq_len(
                                        length(combs_ntaken) / 2)]
                                  }
                                  return(combs_ntaken)
                                  
                                }))
             comb <- sapply(lev_cmb, function(lc) cov %in% lc)
             
             if (coeff_split_type == 'traditional') {
               
               splitpoint <- .select_traditional(values = lev_cmb,
                                                 comb_logical = comb,
                                                 resp = resp)
               
             } else if (coeff_split_type == 'test') {
               
               stat_pval <- apply(comb, 2,
                                  function(q) indep_test(x = q,
                                                         y_dist = resp_dist,
                                                         R = R))
               splitpoint <- .select_splitpoint(values = lev_cmb,
                                                statistic_pvalue = stat_pval)
               
             }
             
           }
           
           # Label levels with 1 if they are in splitpoint, 2 otherwise
           # (and with NA if they do not occur)
           #needed in growtree to split observations using their level
           splitindex <- !(levels(cov) %in% splitpoint)
           splitindex[!(levels(cov) %in% lev)] <- NA_integer_
           splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L
           
         },
         
         fdata      = {
           
           if (split_type == 'coeff') {
             
             comp_idx <- seq_len(dim(new_cov)[2])
             stat_pval <- sapply(comp_idx,
                                 function(i) indep_test(x = new_cov[, i],
                                                        y_dist = resp_dist,
                                                        R = R))
             bselect <- .select_component(statistic_pvalue = stat_pval)
             
             sel_comp <- new_cov[, bselect]
             s  <- sort(sel_comp)
             comb <- sapply(s[-length(s)], function(j) sel_comp <= j)
             
             if (coeff_split_type == 'traditional') {
               
               splitpoint <- .select_traditional(values = s,
                                                 comb_logical = comb,
                                                 resp = resp)
               
             } else if (coeff_split_type == 'test') {
               
               stat_pval <- apply(comb, 2,
                                  function(q) indep_test(x = q,
                                                         y_dist = resp_dist,
                                                         R = R))
               splitpoint <- .select_splitpoint(values = s,
                                                statistic_pvalue = stat_pval)
               
             }
             
           } else if (split_type == 'cluster') {
             
             clustering_objs <- .select_clustering(cov = cov,
                                                   new_cov = new_cov,
                                                   cov_dist = cov_dist)
             splitindex <- clustering_objs$splitindex
             centroids <- clustering_objs$centroids
             
           }
           
         },
         
         list = if (all(sapply(cov, function(obj) attributes(obj)$names)
                        == 'diagram')) {
           
           clustering_objs <- .select_clustering(cov = cov,
                                                 new_cov = new_cov,
                                                 cov_dist = cov_dist)
           splitindex <- clustering_objs$splitindex
           centroids <- clustering_objs$centroids
           
         } else if (all(sapply(cov, class) == 'igraph')) {
           
           if (split_type == 'coeff') {
             
             # Drop non-informative (i.e. all-equal) columns
             new_cov <- new_cov[, !as.logical(apply(new_cov, 2, .zero_range))]
             
             # Control if the df is now void; if so, return 'void'
             if (dim(new_cov)[2] == 0) {
               return(list('void' = TRUE))
             }
             
             comp_idx <- seq_len(dim(new_cov)[2])
             stat_pval <- sapply(comp_idx,
                                 function(i) indep_test(x = new_cov[, i],
                                                        y_dist = resp_dist,
                                                        R = R))
             bselect <- .select_component(statistic_pvalue = stat_pval)
             
             # graph_shell may drop columns, so switch from index to name
             bselect <- as.integer(names(new_cov)[bselect])
             
             sel_comp <- new_cov[[as.character(bselect)]]
             s  <- sort(sel_comp)
             comb <- sapply(s[-length(s)], function(j) sel_comp <= j)
             
             # Check if all columns of comb are equal
             if (all(apply(comb, 2, identical, comb[, 1]))) {
               # if TRUE, the omitted column [position length(s)] is different;
               # when this is the case, set last included column as splitpoint,
               # as it means that breaks is set before last column
               splitpoint <- s[(length(s) - 1)]
             } else if (coeff_split_type == 'traditional') {
               
               splitpoint <- .select_traditional(values = s,
                                                 comb_logical = comb,
                                                 resp = resp)
               
             } else if (coeff_split_type == 'test') {
               
               stat_pval <- apply(comb, 2,
                                  function(q) indep_test(x = q,
                                                         y_dist = resp_dist,
                                                         R = R))
               splitpoint <- .select_splitpoint(values = s,
                                                statistic_pvalue = stat_pval)
               
             }
             
           } else if (split_type == 'cluster') {
             
             clustering_objs <- .select_clustering(cov = cov,
                                                   new_cov = new_cov,
                                                   cov_dist = cov_dist)
             splitindex <- clustering_objs$splitindex
             centroids <- clustering_objs$centroids

           }
         }
  )
  
  split_out <- list()
  if (exists('splitpoint')) split_out$splitpoint <- splitpoint
  if (exists('splitindex')) split_out$splitindex <- splitindex
  if (exists('bselect')) split_out$bselect <- bselect
  if (exists('centroids')) split_out$centroids <- centroids
  return(split_out)
  
}


.select_splitpoint <- function(values, statistic_pvalue) {
  
  stopifnot(identical(rownames(statistic_pvalue), c('Statistic', 'Pvalue')))
  
  if (all(is.na(statistic_pvalue['Pvalue', ])) ||
      length(which(statistic_pvalue['Pvalue', ] ==
                   min(statistic_pvalue['Pvalue', ], na.rm = TRUE))) > 1) {
    
    splitpoint <- values[which.max(statistic_pvalue['Statistic', ])]
    
  } else {
    
    splitpoint <- values[which.min(statistic_pvalue['Pvalue', ])]
    
  }
  
  return(splitpoint)
  
}


.select_component <- function(statistic_pvalue) {
  
  stopifnot(identical(rownames(statistic_pvalue), c('Statistic', 'Pvalue')))
  
  if (all(is.na(statistic_pvalue['Pvalue', ])) ||
      length(which(statistic_pvalue['Pvalue', ] ==
                   min(statistic_pvalue['Pvalue', ], na.rm = TRUE))) > 1) {
    
    bselect <- as.integer(which.max(statistic_pvalue['Statistic', ]))
    
  } else {
    
    bselect <- as.integer(which.min(statistic_pvalue['Pvalue', ]))
    
  }
  
  return(bselect)
  
}
.select_variable <- .select_component


.select_clustering <- function(cov, new_cov, cov_dist) {
  
  stopifnot(identical(typeof(cov), 'list'))
  stopifnot(inherits(cov_dist, 'dist'))
  
  # Initialize splitindex
  N_obs <- length(levels(new_cov))
  splitindex <- rep(NA, N_obs)
  
  # Only two observations: split is trivial
  if (length(cov) == 2) {
    
    splitindex[new_cov] <- c(1, 2)
    
    centroids <- if (class(cov) == 'fdata') {
      list(c1 = cov[1, ], c2 = cov[2, ])
    } else {
      list(c1 = cov[[1]], c2 = cov[[2]])
    }
    
  } else {  # More than two observations: split via PAM

    pam_obj <- cluster::pam(cov_dist, k = 2, diss = TRUE)
    cl_index <- pam_obj$clustering
    splitindex[as.integer(droplevels(new_cov))] <- cl_index
    #as.integer(names(cl_index))?

    medindex1 <- pam_obj$id.med[1]
    c1 <- if (class(cov) == 'fdata') cov[medindex1, ] else cov[[medindex1]]
    medindex2 <- pam_obj$id.med[2]
    c2 <- if (class(cov) == 'fdata') cov[medindex2, ] else cov[[medindex2]]
    centroids <- list(c1 = c1, c2 = c2)
    
  }
  
  return(list(splitindex = splitindex, centroids = centroids))
  
}


.select_traditional <- function(values, comb_logical, resp) {
  
  if (class(resp) == 'factor') {
    
    total_gini <- apply(comb_logical, 2,
                        function(c) {
                          
                          y1 <- resp[c]
                          y2 <- resp[!c]
                          # Gini (alt. def.): G(B_r) = 1 - \sum_k p_{kr}^2
                          g1 <- 1 - sum(prop.table(table(y1)) ^ 2)
                          g2 <- 1 - sum(prop.table(table(y2)) ^ 2)
                          
                          total_gini_c <- g1 + g2
                          return(total_gini_c)
                          
                        })
    
    splitpoint <- values[which.min(total_gini)]
    
  } else {
    
    total_rss <- apply(comb_logical, 2,
                       function(c) {
                         
                         y1 <- resp[c]
                         y2 <- resp[!c]
                         v1 <- if (length(y1) == 1) 0 else var(y1)
                         v2 <- if (length(y2) == 1) 0 else var(y2)
                         n1 <- length(y1)
                         n2 <- length(y2)
                         
                         total_rss_c <- (n1 - 1) * v1 + (n2 - 1) * v2
                         return(total_rss_c)
                         
                       })
    
    splitpoint <- values[which.min(total_rss)]
    
  }
  
  return(splitpoint)
  
}