indep_test <- function(x,
                       y,
                       x_dist = NULL,
                       y_dist = NULL,
                       R = 1000) {
  
  # If distances are not provided, compute them
  if (is.null(x_dist)) x_dist <- dist_comp(x)
  if (is.null(y_dist)) y_dist <- dist_comp(y)
  
  # Distance correlation test
  dct <- energy::dcor.test(x_dist, y_dist, R = R)
  
  # Retrieve and return test statistic and p-value
  stat_pval <- if (!is.na(dct$statistic)) {
    c(dct$statistic, dct$p.value)
  } else{
    c(NA, NA)
  }
  names(stat_pval) <- c('Statistic', 'Pvalue')
  return(stat_pval)
  
}