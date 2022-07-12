#' Distances
#'
#' Compute pairwise distances starting from single objects containing the
#' original univariate observations.
#'
#' @param x Object containing the original univariate observations. Currently available
#'   types and the form they need to have to be correctly recognized are the
#'   following:
#' \itemize{
#'   \item Logical: logical vectors;
#'   \item Numeric: numeric or integer vectors;
#'   \item Nominal: factors;
#'   \item Functions: objects of class \code{"fdata"};
#'   \item Graphs: (lists of) objects of class \code{"igraph"};
#'   \item Persistence diagrams: (lists of) objects with
#'   \code{attributes(x)$names == "diagram"}.
#' } 
#' See Details to find out which distance is used in each case.
#' @param lp Integer specifying which norm should be used to compute the
#'   distances for functional data.
#'   
#'   
#' @details 
#' 
#' The distances used in each case are the following:
#' \itemize{
#'   \item Logical: Euclidean distance, implemented via \code{\link[stats:dist]{dist()}};
#'   \item Numeric: Euclidean distance, implemented via \code{\link[stats:dist]{dist()}};
#'   \item Nominal: Gower's distance, implemented via \code{\link[cluster:daisy]{daisy()}};
#'   \item Functions: \eqn{L^{p}}{L^p}-norm, implemented via
#'   \code{\link[fda.usc:metric.lp]{metric.lp()}} with default options;
#'   \item Graphs: Edge Difference distance (Hammond et al., 2013), implemented via
#'   \code{\link[NetworkDistance:nd.edd]{nd.edd()}};
#'   \item Persistence diagrams: Wasserstein distance, implemented via
#'   \code{\link[TDA:wasserstein]{wasserstein()}} with default options;
#' } 
#'    
#' @returns
#' Object of class \code{"dist"} containing the pairwise distances.
#' 
#' @references 
#' 
#' D. K. Hammond, Y. Gur, and C. R. Johnson (2013). Graph diffusion distance: A
#' difference measure for weighted graphs based on the graph laplacian
#' exponential kernel. In \emph{2013 IEEE Global Conference on Signal and
#' Information Processing}, pages 419-422.
#'
#' @examples
#' 
#' # Number of observations
#' nobs <- 10
#' 
#' ## Logical 
#' obj <- as.logical(rbinom(nobs, 1, 0.5))
#' d <- dist_comp(obj)
#' 
#' ## Integer
#' obj <- rpois(nobs, 5)
#' d <- dist_comp(obj)
#' 
#' ## Numeric
#' obj <- rnorm(nobs)
#' d <- dist_comp(obj)
#' 
#' ## Factors
#' obj <- factor(letters[1:nobs])
#' d <- dist_comp(obj)
#' 
#' ## Functional data
#' obj <- fda.usc::rproc2fdata(nobs, seq(0, 1, len = 100), sigma = 1)
#' d <- dist_comp(obj)
#' 
#' ## Graphs
#' obj <- lapply(1:nobs, function(j) igraph::sample_gnp(100, 0.2))
#' d <- dist_comp(obj)
#' 
#' ## Persistence diagrams
#' x <- lapply(rep(100, nobs), function(np) TDA::circleUnif(np))
#' obj <- lapply(x, TDA::ripsDiag, maxdimension = 1, maxscale = 3)
#' d <- dist_comp(obj)
#' 
#' @export


dist_comp <- function(x,
                      lp = 2) {
  
  # Compute the distance/dissimilarity objects
  dist_obj <- switch(class(x),
                     
                     logical    = dist(x),
                     #needed for split point search when split_type = 'coeff'
                     
                     integer    = dist(x),
                     #objects of class integer are not of class numeric
                     
                     numeric    = dist(x),
                     
                     factor     = cluster::daisy(as.data.frame(x)),
                     
                     fdata      = as.dist(fda.usc::metric.lp(x, lp = lp)),
                     
                     list       = {
                       
                       # Graphs
                       if (all(sapply(x, class) == 'igraph')) {
                         
                         # Unweighted
                         if (all(sapply(x, function(i)
                           is.null(igraph::edge.attributes(i)$weight)))) {
                           adj_data <- lapply(x, igraph::as_adjacency_matrix)
                         } else {  # Weighted
                           adj_data <- lapply(x, function(i) {
                             igraph::as_adjacency_matrix(i, attr = 'weight')
                           })
                         }
                         
                         #d is obtained in the same way in the two cases:
                         d <- NetworkDistance::nd.edd(adj_data)
                         return(d$D)
                         
                         # Persistence diagrams
                       } else if (all(sapply(x, function(x) attributes(x)$names)
                                      == 'diagram')) {
                         wass_fun <- function(i, j) TDA::wasserstein(x[[i]]$diagram,
                                                                     x[[j]]$diagram)
                         vec_wass_fun <- Vectorize(wass_fun)
                         d_idx <- seq_along(x)
                         d <- outer(d_idx, d_idx, vec_wass_fun)
                         return(as.dist(d))
                       }
                       
                     }
  )
  
  # Set observations' names (needed for eforest())
  dist_obj <- usedist::dist_setNames(dist_obj, seq_along(x))

  # Return
  return(dist_obj)
  
}



dist_comp_cl <- function(centroid,
                         x,
                         lp = 2) {
  
  switch(class(x),
         
         fdata  = fda.usc::metric.lp(fdata1 = x, fdata2 = centroid, lp = lp),
         
         list   = {
           
           # Graphs
           if (all(sapply(x, class) == 'igraph')) {
             
             # Unweighted
             if (all(sapply(x, function(i)
               is.null(igraph::edge.attributes(i)$weight)))) {
               adj_data <- lapply(x, igraph::as_adjacency_matrix)
               adj_centroid <- igraph::as_adjacency_matrix(centroid)
             } else {  # Weighted
               adj_data <- lapply(x, function(i) {
                 igraph::as_adjacency_matrix(i, attr = 'weight')
               })
               adj_centroid <- igraph::as_adjacency_matrix(centroid,
                                                           attr = 'weight')
             }
             
             #dist_centroid is obtained in the same way in the two cases:
             dist_centroid <- sapply(adj_data, function(i) {
               d <- NetworkDistance::nd.edd(list(i, adj_centroid))
               return(d$D)
               #since the distances are computed pairwise, d$D is numeric (i.e.,
               #not a 'dist' object anymore)
             })
             return(dist_centroid)
             
             # Persistence diagrams
           } else if (all(sapply(x, function(x) attributes(x)$names)
                          == 'diagram')) {
             wass_fun <- function(x, centroid)
               TDA::wasserstein(x$diagram, centroid$diagram)
             vec_wass_fun <- Vectorize(wass_fun, vectorize.args = 'x')
             return(vec_wass_fun(x, centroid))
           }
           
         }
  )
}
