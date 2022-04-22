growtree <- function(id = 1L,
                     response,
                     covariates,
                     weights,
                     minbucket,
                     alpha,
                     R,
                     split_type,
                     coeff_split_type,
                     p_adjust_method,
                     random_covs) {
  
  # If the observations are less than <2*minbucket>, stop here (otherwise, there
  #would not be enough obs to have at least <minbucket> obs in each kid node)
  if (sum(weights) < 2 * minbucket) {
    return(partynode(id = id))
  }
  
  # Find the best split (variable selection & split point search)
  split <- findsplit(response = response,
                     covariates = covariates,
                     alpha = alpha,
                     R = R,
                     split_type = split_type,
                     coeff_split_type = coeff_split_type,
                     p_adjust_method = p_adjust_method,
                     random_covs = random_covs)
  
  # If no split is found, stop here
  if (is.null(split)) {
    return(partynode(id = id))
  }
  
  # Selected variable index and possibly selected basis index
  varid <- split$varid
  if (!is.null(split$basid)) {
    basid <- split$basid
  }
  
  breaks <- split$breaks
  index <- split$index
  
  # Assign the ids to the observations
  kidids <- c()
  switch(class(covariates$cov[[varid]]),
         
         integer = {
           
           kidids[(which(covariates$cov[[varid]] <= breaks))] <- 1
           kidids[(which(covariates$cov[[varid]] > breaks))] <- 2
           
         },
         
         numeric = {
           
           kidids[(which(covariates$cov[[varid]] <= breaks))] <- 1
           kidids[(which(covariates$cov[[varid]] > breaks))] <- 2
           
         },
         
         factor = {
           
           kidids <- index[covariates$cov[[varid]]]
           #replicate 1 in index for each level in splitpoint; 2 otherwise
           #no need for na.exclude()
           
         },
         
         fdata = {
           
           if (split_type == 'coeff') {
             
             kidids[which(covariates$newcov[[varid]][, basid] <= breaks)] <- 1
             kidids[which(covariates$newcov[[varid]][, basid] > breaks)] <- 2
             
           } else if (split_type == 'cluster') {
             
             kidids <- na.exclude(index)
             
           }
         },
         
         list = if (all(sapply(covariates$cov[[varid]],
                               function(x) attributes(x)$names) == 'diagram')) {
           
           kidids <- na.exclude(index)
           
         } else if (all(sapply(covariates$cov[[varid]], class) == 'igraph')) {
           
           if (split_type == 'coeff') {
             
             kidids[which(covariates$newcov[[varid]][as.character(basid)]
                          <= breaks)] <- 1
             kidids[which(covariates$newcov[[varid]][as.character(basid)]
                          > breaks)] <- 2
             #recall that basid is the column's name, so [as.character(basid)]
             
           } else if (split_type == 'cluster') {
             
             kidids <- na.exclude(index)
             
           }
         }
  )
  
  # Initialization of the kid nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  
  # Give birth to the kid nodes
  for (kidid in seq_along(kids)) {
    
    # Select observations for the current kid node
    w <- weights
    w[kidids != kidid] <- 0
    
    # For less than <minbucket> observations, stop here
    if (sum(w) < minbucket) {
      return(partynode(id = id))
    }
    
    # Previous maximum id (to later set the id for the current kid node)
    if (kidid > 1) {
      prev_id <- max(nodeids(kids[[kidid - 1]]))
    } else {
      prev_id <- id
    }
    
    # Recursion on this kid node: update covariates and response
    covariates_updated <- list()
    covariates_updated$cov <- lapply(covariates$cov,
                                     function(cov) subset(cov, as.logical(w)))
    covariates_updated$newcov <- lapply(covariates$newcov, function(newcov)
      subset(newcov, as.logical(w)))
    covariates_updated$dist <- lapply(covariates$dist, function(cov_dist)
      usedist::dist_subset(cov_dist, which(w == 1)))
    response_updated <- list()
    response_updated$response <- subset(response$response, as.logical(w))
    response_updated$response_dist <-
      usedist::dist_subset(response$response_dist, which(w == 1))
    
    # Actual recursion
    kids[[kidid]] <- growtree(id = as.integer(prev_id + 1),
                              response = response_updated,
                              covariates = covariates_updated,
                              weights = rep(1L, sum(w, na.rm = TRUE)),
                              minbucket = minbucket,
                              alpha = alpha,
                              R = R,
                              split_type = split_type,
                              coeff_split_type = coeff_split_type,
                              p_adjust_method = p_adjust_method,
                              random_covs = random_covs)
  }
  
  # Return the nodes (i.e. the split rules)
  return(partynode(id = as.integer(id),
                   split = split,
                   kids = kids,
                   info = list(pvalue = min(info_split(split)$pvalue,
                                            na.rm = TRUE))
  ))
}