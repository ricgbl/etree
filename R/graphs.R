graph_shell <- function(graph_list,
                        max_shell = NULL,
                        predicting = FALSE) {
  
  # Number of observations (graphs)
  n_graphs <- length(graph_list)
  
  # Shell distribution for each graph
  table_shell <- lapply(graph_list, function(g) table(brainGraph::s_core(g)))
  
  # Observed maximum shell index
  obs_max_shell <- do.call(max, lapply(table_shell,
                                       function(s) as.integer(names(s))))
  
  # If max_shell is given, go for it; otherwise, set obs_max_shell as max_shell
  if (is.null(max_shell)) max_shell <- obs_max_shell
  
  # Column names for the shell df
  col_names <- as.character(seq(0, max_shell, 1))
  
  # Shell df inizialization
  all_shell_df <- data.frame(matrix(data = 0L,
                                    nrow = n_graphs,
                                    ncol = length(col_names)))
  colnames(all_shell_df) <- col_names
  
  # Fill in with the actual shell distibutions
  invisible(sapply(seq_len(n_graphs), function(i) {
    shells <- names(table_shell[[i]])
    cols <- intersect(col_names, shells)
    all_shell_df[i, cols] <<- table_shell[[i]][cols]
  }))
  # better a for cycle?
  # for(i in 1:n.graphs){
  #   cols <- names(table.shell[[i]])
  #   all.shell.df[i, cols] = table.shell[[i]][cols]
  # }
  
  # Ignore non-informative columns only when not in predict
  if (isFALSE(predicting)) {
    # Update all_shell_df ignoring non-informative columns
    all_shell_df <- all_shell_df[, !as.logical(apply(all_shell_df,
                                                     2,
                                                     .zero_range))]
  }
  
  # Return the final shell df
  return(all_shell_df)
  
}


# Function to test if all elements of a given vector are equal for tol provided
.zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  # if only one obs, equality cannot be tested -> return FALSE
  if (length(x) == 1) {
    return(FALSE)
  }
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}