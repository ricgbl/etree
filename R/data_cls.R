#' Classification toy dataset
#' 
#' A simple dataset containing simulated values for a nominal response variable
#' and four covariates of both mixed and partially structured type. The data
#' generation process is based on Example 4.7 (''Signal shape classification'',
#' pages 73-77) from Saito (1994).
#' 
#' @format List with two elements: \code{covs}, which is a list containing the
#'   covariates, and \code{resp}, which is a factor of length 150 representing
#'   the response variable. The response variable is divided into three classes
#'   whose labels are cylinder (\code{Cyl}), bell (\code{Bel}) and funnel
#'   (\code{Fun}). The four covariates in \code{covs} all have length 150 and
#'   are characterized as follows:
#' \itemize{
#' \item Nominal: \code{Cyl} observations are given level 1 with probability 0.8
#' and levels 2 and 3 with probability 0.1 each, \code{Bel} observations are
#' given level 2 with probability 0.8 and levels 1 and 3 with probability 0.1
#' each, \code{Fun} observations are given level 3 with probability 0.8 and
#' levels 1 and 2 with probability 0.1 each;
#' \item Numeric: coefficients for one of the basis used to perform the
#' B-splines expansion of the curves that are in turn specified as in Saito
#' (1994);
#' \item Functional: curves as specified in Saito (1994);
#' \item Graphs: Erd\"{o}s-R\'{e}nyi graphs with connection probability 0.10 for
#' \code{Cyl} observations, 0.125 for \code{Bel} observations, 0.15 for
#' \code{Fun} observations.
#' }
#' 
#' @references
#'  
#' Saito, N. (1994). Local feature extraction and its applications using a
#' library of bases (Doctoral dissertation, Yale University).
#' 
"data_cls"