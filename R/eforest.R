#' Energy Forests
#'
#' Fits an Energy Forest, in the form of either a bagging of Energy Trees or a
#' Random Energy Forest, depending on the value of the \code{random_covs}
#' parameter.
#'
#' @param ntrees Number of Energy Trees to be built, i.e., the number of
#'   bootstrap samples to be generated and used for fitting.
#' @param ncores Number of cores to use, i.e., at most how many child processes
#'   will be run simultaneously. Must be exactly 1 on Windows (which uses the
#'   master process). \code{ncores} corresponds to \code{mc.cores} in
#'   \code{\link[parallel:mclapply]{mclapply()}}, which is actually used to grow the single
#'   Energy Trees in a parallel fashion.
#' @param perf_metric Performance metric that is used to compute the Out-Of-Bag
#'   score. If \code{NULL}, default choices are used: Balanced Accuracy for
#'   classification and Root Mean Square Percentage Error for regression. See
#'   Details for further information and possible alternatives.
#' @param random_covs Size of the random subset of covariates to choose from at
#'   each split. If set to \code{NULL} (default), all the covariates are
#'   considered each time, resulting in a bagging of Energy Trees. When
#'   \code{random_covs} is an integer greater than 1 and less than the total
#'   number of covariates, the model is a Random Energy Forest.
#' @param verbose Logical indicating whether to print a one-line notification
#'   for the conclusion of each tree's fitting process.
#' @inheritParams etree
#'
#' @details
#' 
#' \code{eforest()} generates \code{ntrees} bootstrap samples and then calls
#' \code{\link[etree:etree]{etree()}} on each of them. Then, it computes the Out-Of-Bag (OOB)
#' score using the performance metric defined through \code{perf_metric}.
#' 
#' For classification, possible values of \code{perf_metric} are \code{"BAcc"}
#' and \code{"WBAcc"}. Both are general enough to be used in multiclass
#' classification problems, still producing sensible results in the case of
#' binary classification. The two options are based on the calculation of a
#' ground performance metric, the Balanced Accuracy, which is defined as the
#' arithmetic mean between Sensitivity and Specificity. In this framework,
#' Balanced Accuracy is computed using a "One vs. All" approach, i.e.,
#' considering one class at a time: positive instances are those belonging to
#' that class, and negatives are the ones belonging to any other class. Then,
#' the "One vs. All" Balanced Accuracies obtained by considering each class must
#' be averaged. When \code{perf_metric = "BAcc"} (default for classification
#' tasks), the average is arithmetic. When \code{perf_metric = "WBAcc"}, the
#' average is weighted using class sizes, hence giving more importance to the
#' "One vs. All" Balanced Accuracy of larger classes.
#' 
#' For regression, the default value of \code{perf_metric} is \code{"RMSPE"},
#' namely, Root Mean Square Percentage Error. Other available options are
#' \code{c("MAE", "MAPE", "MedianAE", "MedianAPE", "MSE", "NRMSE", "RAE",
#' "RMSE", "RMLSE")}. Each of these name points to the corresponding homonym
#' function from the package \code{\link[MLmetrics]{MLmetrics}}, whose
#' documentation provides more information about their definition.
#'
#' @section Value:
#' Object of class \code{"eforest"} with three elements: 1) \code{ensemble},
#' which is itself a list gathering all the fitted trees; 2) \code{oob_score},
#' an object of class \code{"numeric"} representing the OOB score computed using
#' the performance metric defined through \code{perf_metric}; 3)
#' \code{perf_metric}, an object of class \code{"character"} returning the
#' performance metric used for computations.
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
#' ## Response variable(s)
#' resp_reg <- cov_num ^ 2
#' y <- round((cov_num - min(cov_num)) / (max(cov_num) - min(cov_num)), 0)
#' resp_cls <- factor(y)
#' 
#' ## Regression ##
#' eforest_fit <- eforest(response = resp_reg, covariates = cov_list, ntrees = 20)
#' print(eforest_fit$ensemble[[1]])
#' plot(eforest_fit$ensemble[[1]])
#' mean((resp_reg - predict(eforest_fit)) ^ 2)
#' 
#' ## Classification ##
#' eforest_fit <- eforest(response = resp_cls, covariates = cov_list, ntrees = 20)
#' print(eforest_fit$ensemble[[20]])
#' plot(eforest_fit[[20]])
#' table(resp_cls, predict(eforest_fit))
#' }
#'
#' @export

eforest <- function(response,
                    covariates,
                    weights = NULL,
                    ntrees = 100,
                    ncores = 1L,
                    minbucket = 1,
                    alpha = 1,
                    R = 500,
                    split_type = 'cluster',
                    coeff_split_type = 'test',
                    p_adjust_method = 'fdr',
                    perf_metric = NULL,
                    random_covs = 'auto',
                    verbose = FALSE) {
  
  # Check whether covariates is a list
  if (!is.list(covariates)) {
    stop("Argument 'covariates' must be provided as a list")
  }
  
  # Check that response is factor or numeric
  if (!is.factor(response) &&
      !is.numeric(response)) {
    stop("Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  }
  
  # Check that minbucket is positive and ensure it is integer
  stopifnot((minbucket > 0))
  minbucket <- as.integer(minbucket)
  
  # Check that alpha is between 0 and 1
  stopifnot((alpha >= 0) && (alpha <= 1))
  
  # Check that ntrees is positive and ensure it is integer
  stopifnot((ntrees > 0))
  ntrees <- as.integer(ntrees)
  
  # If perf_metric is NULL, set it to default choices
  if (is.null(perf_metric)) {
    if (is.factor(response)) perf_metric <- 'BAcc' else
      if (is.numeric(response)) perf_metric <- 'RMSPE'
  } else {
    # If perf_metric is given, check that it has an acceptable value
    if (is.factor(response)) {
      if (isFALSE(perf_metric %in% c('BAcc', 'WBAcc'))) {
        stop("The value provided for 'perf_metric' is not available for this
             task")
      }
    } else if (is.numeric(response)) {
      if (isFALSE(perf_metric %in% c('MAPE', 'RMSPE', 'NRMSE', 'MAE',
                                     'MedianAE', 'MedianAPE', 'MSE', 'RAE',
                                     'RMSE', 'RMLSE'))) {
        stop("The value provided for 'perf_metric' is not available for this
             task")
      }
    }
  }
  
  # If random_covs is 'auto', set it to traditional values
  if (random_covs == 'auto') {
    if (is.factor(response)) {
      random_covs <- floor(sqrt(length(covariates)))
    } else if (is.numeric(response)) {
      random_covs <- floor(length(covariates) / 3)
    } 
  }
  
  # Number of observations
  nobs <- length(response)
  
  # New list of covariates
  newcovariates <- .create_newcov(covariates = covariates,
                                  response = response,
                                  split_type = split_type)
  
  # Distances
  cov_distance <- lapply(covariates, dist_comp)
  
  # Large list with covariates, newcovariates and distances
  covariates_large <- list('cov' = covariates,
                           'newcov' = newcovariates,
                           'dist' = cov_distance)
  
  # Large list with response and the corresponding distances
  response_large <- list('response' = response,
                         'response_dist' = dist_comp(response))
  
  # Generate B bootstrap samples
  set.seed(12345)
  boot_idx <- lapply(1:ntrees,
                     function(b) sample.int(nobs, replace = TRUE))
  
  # Energy Tree fits for each bootstrap sample
  etree_boot_fits <- parallel::mclapply(boot_idx, function(b_i) {
    
    # Covariates and response for this bootstrap sample
    boot_cov_large <- lapply(covariates_large[1:2],
                             function(cl) lapply(cl, function(cov) {
                               if (class(cov) == 'data.frame') cov[b_i, ]
                               else cov[b_i]
                             }
                             ))
    # Re-index newcov only if using 'cluster'
    if (split_type == 'cluster') {
      boot_cov_large$newcov[[1]] <- factor(1:nobs)
      boot_cov_large$newcov[[2]] <- factor(1:nobs)
    }
    boot_cov_large$dist <- lapply(covariates_large[[3]],
                                  function(cov_dist) {
                                    boot_dist <- usedist::dist_subset(cov_dist,
                                                                      b_i)
                                    boot_dist <-
                                      usedist::dist_setNames(boot_dist, 1:nobs)
                                    return(boot_dist)
                                  })
    resp_dist <- usedist::dist_subset(response_large$response_dist, b_i)
    resp_dist <- usedist::dist_setNames(resp_dist, 1:nobs)
    boot_resp_large <- list('response' = response[b_i],
                            'response_dist' = resp_dist)
    
    # Energy Trees fit
    e_fit <- etree(response = boot_resp_large,
                   covariates = boot_cov_large,
                   weights = weights,
                   minbucket = minbucket,
                   alpha = alpha,
                   R = R,
                   split_type = split_type,
                   coeff_split_type = coeff_split_type,
                   p_adjust_method = p_adjust_method,
                   random_covs = random_covs)
    
    # Remove .Environment attribute from terms
    attr(e_fit$terms, '.Environment') <- NULL
    
    # Print counter if verbose is TRUE
    if (isTRUE(verbose)) {
      if (isTRUE(exists('counter'))) counter <<- counter + 1 else counter <<- 1
      print(paste('Fitted Energy Tree n.', counter))
    }
    
    # Return
    return(e_fit)
    
  },
  mc.cores = ncores,
  mc.set.seed = TRUE)
  
  # Delete counter if verbose is TRUE
  if (isTRUE(verbose)) rm(counter, pos = 1)
  
  # For each obs, predicted response with its own OOB trees
  oob_pred_resp <- lapply(1:nobs, function(i) {
    
    # Covariates for observation i (to be used as newdata in predict())
    covs_i <- lapply(covariates, function(cov) cov[i])
    
    # Indices of bootstrap samples for which the observation is OOB
    is_oob <- sapply(boot_idx, function(b) !(i %in% b))
    is_oob_idx <- which(is_oob)
    
    # Predict response for this obs only with trees with indices in is_oob_idx
    oob_pred_resp <- sapply(is_oob_idx, function(o) {
      predict(etree_boot_fits[[o]],
              newdata = covs_i)
    })
    
    # Return predicted response
    return(oob_pred_resp)
    
  })
  
  # Predicted responses and OOB performance metric calculation
  if (is.factor(response)) {
    
    ## Classification ##
    
    # Predicted response: majority voting rule
    pred_resp <- factor(sapply(oob_pred_resp,
                               function(i) names(which.max(table(i)))
    ))
    
    # OOB performance metric (measured via BAcc or WBAcc)
    if (perf_metric == 'BAcc' || perf_metric == 'WBAcc') {
      
      # Balanced Accuracy for each class (each given by (TP/P + TN/N) / 2)
      bal_accs <- .comp_bal_accs(y_pred = pred_resp,
                                 y_true = response)
      
      # OOB score
      #Balanced Accuracy -> arithmetic mean
      #Weighted Balanced Accuracy -> weighted mean (with class sizes as weights)
      oob_score <- switch(perf_metric,
                          BAcc = mean(bal_accs),
                          WBAcc = sum(bal_accs * table(response)) /
                            length(response))
      
    }
    
  } else if (is.numeric(response)) {
    
    ## Regression ##
    
    # Predicted response: average
    pred_resp <- sapply(oob_pred_resp, mean)
    
    # OOB performance metric: various choices (default is 'RMSPE')
    oob_score <-
      switch(perf_metric,
             MAPE = MLmetrics::MAPE(pred_resp, response),
             RMSPE = MLmetrics::RMSPE(pred_resp, response),
             NRMSE = MLmetrics::RMSE(pred_resp, response) / mean(response),
             MAE = MLmetrics::MAE(pred_resp, response),
             MedianAE = MLmetrics::MedianAE(pred_resp, response),
             MedianAPE = MLmetrics::MedianAPE(pred_resp, response),
             MSE = MLmetrics::MSE(pred_resp, response),
             RAE = MLmetrics::RAE(pred_resp, response),
             RMSE = MLmetrics::RMSE(pred_resp, response),
             RMLSE = MLmetrics::RMSLE(pred_resp, response))
    
  }
  
  # Create object with etree_boot_fits, oob_score, and perf_metric
  eforest_obj <- list(ensemble =  etree_boot_fits,
                      oob_score = oob_score,
                      perf_metric = perf_metric)
  class(eforest_obj) <- 'eforest'
  
  # Return eforest object
  return(eforest_obj)
  
}


#' Predictions for Energy Forests
#'
#' Compute predictions for objects of class \code{"eforest"} (i.e., as returned
#' by \code{\link[etree:eforest]{eforest()}}).
#'
#' @param object A fitted Energy Forest of class \code{"eforest"}.
#' @param newdata Optional set of new covariates used to make predictions. Must
#'   be provided as a list, where each element is a different variable.
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
#' \code{newdata} is omitted, fitted values of individual trees are somehow
#' combined (see Details) and returned.
#' @param ... Additional arguments.
#'
#' @details
#' The \code{predict()} method for \code{"eforest"} objects computes predictions
#' for Energy Forests as returned by \code{\link[etree:eforest]{eforest()}}.
#' Predictions are based either on the fitted values (if \code{newdata} is
#' \code{NULL}) or on the new set of covariates (when \code{newdata} is
#' provided). In both cases, each tree in \code{object$ensemble} is used to make
#' predictions by calling \code{\link[etree:predict.etree]{predict()}} on it
#' (with the same specification of \code{newdata}). Then, individual trees'
#' predictions for any single observation are combined by majority voting rule
#' for classification or by arithmetic mean for regression.
#' 
#' @section Value:
#' Predictions, in the form of a factor for classification or as a numeric
#' vector for regression.
#'
# @examples  
#'
#' @method predict eforest  
#' @export

predict.eforest <- function(object, newdata = NULL, ...) {
  
  # Individual predictions with base learners
  #(newdata check, split retrieval, newcov computation are all done in predict)
  ind_pred_resp <- sapply(object$ensemble, function(tree) {
    predict(tree, newdata = newdata)
  })
  
  # Predict response, differently based on the type of problem (CLS or REG)
  response <- object$ensemble[[1]]$fitted$`(response)`
  if (is.factor(response)) {
    
    # Majority voting rule
    pred_resp <- apply(ind_pred_resp, 1, function(i) names(which.max(table(i))))
    pred_resp <- factor(pred_resp)
    
  } else if (is.numeric(response)) {
    
    # Average
    pred_resp <- apply(ind_pred_resp, 1, mean)
    
  }
  
  # Return predicted response
  return(pred_resp)
  
}


.comp_bal_accs <- function(y_pred,
                          y_true) {
  
  bal_accs <- sapply(levels(y_true),
                     function(lev) {
                       # Sensitivity
                       true_pos <- sum(y_pred == lev &
                                         y_true == lev)
                       pos <- sum(y_true == lev)
                       sens <- true_pos / pos
                       # Specificity
                       true_neg <- sum(y_pred != lev &
                                         y_true != lev)
                       neg <- sum(y_true != lev)
                       spec <- true_neg / neg
                       # Balanced Accuracy (with the current class as pos)
                       bal_acc_lev <- (sens + spec) / 2
                       # Return Balanced Accuracy
                       return(bal_acc_lev)
                     })
  
  return(bal_accs)
  
  
  
}
