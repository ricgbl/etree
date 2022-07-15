#' Energy Tree
#'
#' Fits an Energy Tree for classification or regression.
#'
#' @param response Response variable, an object of class either
#'   \code{"factor"} or \code{"numeric"} (for classification and regression,
#'   respectively).
#' @param covariates Set of covariates. Must be provided as a list, where
#'   each element is a different variable. Currently available types and the
#'   form they need to have to be correctly recognized are the following:
#' \itemize{
#'   \item Numeric: numeric or integer vectors;
#'   \item Nominal: factors;
#'   \item Functions: objects of class \code{"fdata"};
#'   \item Graphs: (lists of) objects of class \code{"igraph"}.
#   \item Persistence diagrams: (lists of) objects with
#   \code{attributes(x)$names == "diagram"}.
#' } 
#' Each element (i.e., variable) in the covariates list must have the same
#' \code{length()}, which corresponds to the sample size.
#' @param weights Optional vector of non-negative integer-valued weights to
#'   be used in the fitting process. If not provided, all observations are
#'   assumed to have weight equal to 1.
#' @param minbucket Positive integer specifying the minimum number of
#'   observations that each terminal node must contain. Default is 5.
#' @param alpha Nominal level controlling the probability of type I error in the
#'   Energy tests of independence used for variable selection. Default is 0.05.
#' @param R Number of replicates employed to approximate the sampling
#'   distribution of the test statistic in every Energy test of independence.
#'   Default is 1000.
#' @param split_type Splitting method used when the selected covariate is
#'   structured. It has two possible values: \code{"coeff"} for feature vector
#'   extraction, and \code{"cluster"} for clustering. See Details for further
#'   information.
#' @param coeff_split_type Method to select the split point for the chosen
#'   component when the selected covariate is structured and \code{split_type =
#'   "coeff"}. It has two possible values: \code{"test"}, in which case Energy
#'   tests of independence are used, and \code{"traditional"}, to employ
#'   traditional methods (Gini index for classification and RSS for regression).
#'   See Details for further information.
#' @param p_adjust_method Multiple-testing adjustment method for P-values,
#'   which can be set to any of the values provided by
#'   \code{\link[stats]{p.adjust.methods}}. Default is \code{"fdr"} for False
#'   Discovery Rate.
# @param supervised whether (\code{TRUE}) or not (\code{FALSE}) to use the
#   supervised version of the model. Default is \code{TRUE}.
#' @param random_covs Size of the random subset of covariates to choose from
#'   at each split. If set to \code{NULL} (default), all the covariates are
#'   considered each time.
#'
#' @details
#' 
#' \code{etree()} is the main function of the homonym package. It allows
#' implementing Energy Trees by simply specifying the response variable, the set
#' of covariates, and possibly some other parameters. The function is specified
#' in the same way regardless of the task type: the choice between
#' classification and regression is automatically made depending on the nature
#' of the response variable.
#' 
#' Energy Trees (Giubilei et al., 2022) are a recursive partitioning tree-based 
#' model built upon
#' Conditional Trees (Hothorn et al., 2006). At each step of Energy Trees'
#' iterative procedure, an Energy test of independence (Szekely et al., 2007) is
#' performed between the response variable and each of the J covariates. If the
#' test of global independence (defined as the intersection of the J tests of
#' partial independence) is not rejected at the significance level set by
#' \code{alpha}, the recursion is stopped; otherwise, the covariate most
#' associated with the response in terms of P-value is selected for splitting.
#' When the covariate is traditional (i.e, numeric or nominal), an Energy test
#' of independence is performed for each possible split point, and the one
#' yielding the strongest association with the response is chosen. When the
#' selected covariate is structured, the split procedure is defined by the value
#' of \code{split_type}, and possibly by that of \code{coeff_split_type}.
#' 
#' \code{split_type} specifies the splitting method for structured covariates.
#' It has two possible values:
#' \itemize{
#' \item \code{"coeff"}: in this case, feature vector extraction is used to
#' transform the structured selected covariate into a set of numeric components
#' using a representation that is specific to its type. Available
#' transformations of such a kind are cubic B-spline expansions for functional
#' data and shell distributions (Carmi et al., 2007) for graphs - obtained
#' through k-cores (Seidman, 1983), s-cores (Eidsaa and Almaas, 2013), and
#' d-cores (Giatsidis et al., 2013), for binary, weighted, and directed graphs,
#' respectively. Then, the component most associated with the response is
#' selected using Energy tests of independence (Szekely et al., 2007), and the
#' split point for that component is chosen using the method defined by
#' \code{coeff_split_type};
#' \item \code{"cluster"}: in this case, the observed values for the structured
#' selected covariate are used within a Partitioning Around Medoids (Kaufmann
#' and Rousseeuw, 1987) step to split observations into the two kid nodes.
#' Medoids calculation and units assignment are performed using
#' \code{\link[cluster:pam]{pam()}}. Distances are specific to each type of
#' variable (see \code{\link[etree:dist_comp]{dist_comp()}} for details).
#' }
#' 
#' \code{coeff_split_type} defines the method to select the split point for the
#' chosen component of the selected structured covariate if and only if
#' \code{split_type = "coeff"}. It has two possible values: 
#' \itemize{ 
#' \item \code{"test"}: an Energy test of independence (Szekely et al., 2007) is
#' performed for each possible split point of the chosen component, and the one
#' yielding the strongest association with the response is selected;
#' \item \code{"traditional"}: the split point for the chosen component is
#' selected as the one minimizing the Gini index (for classification) or the RSS
#' (for regression) in the two kid nodes.
#' }
#'
#' @returns
#' An object of class \code{"etree"}, \code{"constparty"}, and \code{"party"}.
#' It stores all the information about the fitted tree. Its elements can be
#' individually accessed using the \code{$} operator. Their names and content
#' are the following:
#' \itemize{
#' \item \code{node}: a \code{\link[partykit]{partynode}} object representing
#' the basic structure of the tree;
#' \item \code{data}: a \code{list} containing the data used for the fitting
#' process. Traditional covariates are included in their original form, while
#' structured covariates are stored in the form of components if
#' \code{split_type = "coeff"} or as a \code{factor} whose levels go from 1 to 
#' the total number of observations if \code{split_type = "cluster"};
#' \item \code{fitted}: a \code{data.frame} whose number of rows coincides with 
#' the sample size. It includes the fitted terminal node identifiers (in 
#' \code{"(fitted)"}) and the response values of all observations (in 
#' \code{"(response)"});
#' \item \code{terms}: a \code{\link[stats]{terms}} object;
#' \item \code{names} (optional): names of the nodes in the tree. They can be
#' set using a \code{character} vector: if its length is smaller than the number
#' of nodes, the remaining nodes have missing names; if its length is larger,
#' exceeding names are ignored.
#' }
#' 
#' @seealso 
#' \code{\link[partykit:ctree]{ctree()}} for the \code{partykit} implementation of
#' Conditional Trees (Hothorn et al., 2006).
#' 
#' @references 
#' 
#' R. Giubilei, T. Padellini, P. Brutti (2022). Energy Trees: Regression and 
#' Classification With Structured and Mixed-Type Covariates. arXiv preprint.
#' https://arxiv.org/pdf/2207.04430.pdf.
#' 
#' S. Carmi, S. Havlin, S. Kirkpatrick, Y. Shavitt, and E. Shir (2007). A model
#' of internet topology using k-shell decomposition. \emph{Proceedings of the
#' National Academy of Sciences}, 104(27):11150-11154.
#' 
#' M. Eidsaa and E. Almaas (2013). S-core network decomposition: A
#' generalization of k-core analysis to weighted networks. \emph{Physical Review
#' E}, 88(6):062819.
#' 
#' C. Giatsidis, D. M. Thilikos, and M. Vazirgiannis (2013). D-cores: measuring
#' collaboration of directed graphs based on degeneracy. \emph{Knowledge and
#' information systems}, 35(2):311-343.
#' 
#' T. Hothorn, K. Hornik, and A. Zeileis (2006). Unbiased recursive
#' partitioning: A conditional inference framework. \emph{Journal of
#' Computational and Graphical Statistics}, 15(3):651-674.
#' 
#' L. Kaufmann and P. Rousseeuw (1987). Clustering by means of medoids.
#' \emph{Data Analysis based on the L1-Norm and Related Methods}, pages 405-416.
#' 
#' S. B. Seidman (1983). Network structure and minimum degree. \emph{Social
#' networks}, 5(3):269-287.
#' 
#' G. J. Szekely, M. L. Rizzo, and N. K. Bakirov (2007). Measuring and testing
#' dependence by correlation of distances. \emph{The Annals of Statistics},
#' 35(6):2769-2794.
#'
#' @examples
#' \donttest{
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
#' etree_fit <- etree(response = resp_reg, covariates = cov_list)
#' print(etree_fit)
#' plot(etree_fit)
#' mean((resp_reg - predict(etree_fit)) ^ 2)
#' 
#' ## Classification ##
#' etree_fit <- etree(response = resp_cls, covariates = cov_list)
#' print(etree_fit)
#' plot(etree_fit)
#' table(resp_cls, predict(etree_fit))
#' }
#'
#' @export


etree <- function(response,
                  covariates,
                  weights = NULL,
                  minbucket = 5,
                  alpha = 0.05,
                  R = 1000,
                  split_type = "coeff",
                  coeff_split_type = "test",
                  p_adjust_method = "fdr",
#                 supervised = TRUE,
                  random_covs = NULL) {
  
  # Check whether covariates is a list
  if (!is.list(covariates)) {
    stop("Argument 'covariates' must be provided as a list")
  }
  
  # Check that minbucket is positive and ensure it is integer
  stopifnot((minbucket > 0))
  minbucket <- as.integer(minbucket)
  
  # Check that alpha is between 0 and 1
  stopifnot((alpha >= 0) && (alpha <= 1))
  
  # If the case weights are not provided, they are all initialized as 1
  if (is.null(weights)) {
    if (isTRUE(identical(names(response), c('response', 'response_dist')))) {
      weights <- rep(1L, as.numeric(length(response$response)))
    } else {
      weights <- rep(1L, as.numeric(length(response)))
    }
  }
  
  # Large list both for covariates and for response
  #(get them, if already computed in eforest; otherwise, compute them)
  if (isTRUE(identical(names(covariates), c('cov', 'newcov', 'dist')))) {
    
    # Check that response is factor or numeric
    if (!is.factor(response$response) &&
        !is.numeric(response$response)) {
      stop("Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
    }
    
    # Large list with covariates, newcovariates and distances
    covariates_large <- covariates
    
    # Large list with response and the corresponding distances
    response_large <- response
    
    # New covariates and response (needed divided to build the df used by party)
    newcovariates <- covariates_large$newcov
    response <- response_large$response
    
  } else {
    
    # Check that response is factor or numeric
    if (!is.factor(response) &&
        !is.numeric(response)) {
      stop("Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
    }
    
    # New covariates
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
    
  }
  
  # Grow the tree (find the split rules)
  nodes <- growtree(id = 1L,
                    response = response_large,
                    covariates = covariates_large,
                    weights = weights,
                    minbucket = minbucket,
                    alpha = alpha,
                    R = R,
                    split_type = split_type,
                    coeff_split_type = coeff_split_type,
                    p_adjust_method = p_adjust_method,
                    random_covs = random_covs)
  
  # Actually perform the splits
  fitted_obs <- fitted_node(nodes, data = newcovariates)
  
  # Return a rich constparty object
  obj <- party(nodes,
               data = newcovariates,
               fitted = data.frame("(fitted)" = fitted_obs,
                                   "(response)" = response,
                                   #if (supervised) response else fitted_obs,
                                   check.names = FALSE),
               terms = terms(response ~ ., data = newcovariates))
  etree_obj <- as.constparty(obj)
  attr(etree_obj, 'split_type') <- split_type     #used in predict.party
  
  class(etree_obj) <- c('etree', 'constparty', 'party')
  return(etree_obj)
  
}


.create_newcov <- function(covariates, response, split_type) {
  
  newcovariates <- lapply(covariates, function(j) {
    
    switch(class(j),
           
           fdata = {
             
             if (split_type == "coeff") {
               
               fdata_est <- fda.usc::optim.basis(j, numbasis =
                                                   floor(seq(4,
                                                             ncol(j) / 2,
                                                             len = 10)))
               #seq starts from 4 as it is the smallest ok number (norder is 4,
               #and nbasis has to be >= norder -- cf. fda::create.bspline.basis)
               coefs <- fda.usc::fdata2fd(fdata_est$fdata.est,
                                          type.basis = "bspline",
                                          nbasis = fdata_est$numbasis.opt)$coefs
               newcov <- data.frame(t(coefs))
               names(newcov) <- seq_along(names(newcov))
               
             } else if (split_type == "cluster") {
               
               newcov <- as.factor(seq_along(response))
               
             }
             
             attr(newcov, 'cov_type') <- 'fdata'
             return(newcov)
             
           },
           
           list = if (all(sapply(j, function(x) attributes(x)$names)
                          == 'diagram')) {
             
             newcov <- as.factor(seq_along(response))
             attr(newcov, 'cov_type') <- 'diagram'
             return(newcov)
             
           } else if (all(sapply(j, class) == 'igraph')) {
             
             if (split_type == "coeff") {
               
               newcov <- graph_shell(j)
               
             } else if (split_type == "cluster") {
               
               newcov <- as.factor(seq_along(response))
               
             }
             
             attr(newcov, 'cov_type') <- 'graph'
             return(newcov)
             
           },
           
           # In all other cases:
           return(j)
           
    )
  })
  
  # Covariates name
  if (!is.null(names(covariates))) {
    names(newcovariates) <- names(covariates)
    #control if any name is void, i.e. if it is ''
    no_name <- which(sapply(names(covariates), function(n) n == '',
                            USE.NAMES = FALSE))
    names(newcovariates) <- replace(names(newcovariates),
                                    no_name,
                                    as.factor(seq_along(no_name)))
  } else {
    warning('No names available for covariates. Numbers are used instead.')
    names(newcovariates) <- seq_along(newcovariates)
  }
  
  # Return newcovariates
  return(newcovariates)
  
}
