#' @importFrom ripserr vietoris_rips
#' @importFrom simplextree simplex_tree expand
NULL


`%||%` <- getFromNamespace("%||%", "ggplot2")
ensure_nonempty_data <- getFromNamespace("ensure_nonempty_data", "ggplot2")


# Wrapper around ripserr:: functions for calculating persistent homologies
data_to_persistent_homology <- function(df, diameter_max = -1, k = 1, complex = "Rips") {

  # Need to implement other complexes (Cech, Alpha, ...)
  if (complex != "Rips") stop("Only Rips complex is currently implemented")

  d <- as.matrix(df)

  res <- ripserr::vietoris_rips(d, dim = k, threshold = diameter_max, p = 2L, return_format = "df")
  names(res) <- c("dim", "birth", "death")

  res

}


# Returns simplextree representing complex corresponding to diameter w/ simplexes up to dim k
data_to_simplextree <- function(df, diameter, k = 10, complex = "Rips") {

  # Need to implement other complexes (Cech, Alpha, ...)
  if (complex != "Rips") stop("Only Rips complex is currently implemented")

  # Find edges given diameter:
  edges <- t(proximate_pairs(df, diameter))
  edges <- as.data.frame(edges)
  edges <- as.list(edges)

  # Construct the 1-skeleton of the complex as a simplextree
  st <- simplextree::simplex_tree(edges)

  # For the Rips complex, just return the flag complex:
  # w/ simplexes of maximal dim k
  simplextree::expand(st, k = k)

}

proximate_pairs <- function(data, diameter) {

  distances <- as.matrix(stats::dist(data[, c("x", "y")]))
  which(distances < diameter & upper.tri(distances), arr.ind = TRUE)

}
