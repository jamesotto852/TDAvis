#' Persistence homologies from point clouds
#'
#' Compute persistence homologies and plot as barcode charts and persistence diagrams.
#' Accepts either a single point cloud of arbitrary dimension via the `point_cloud` argument
#' or multiple data sets via the `point_data` aesthetic mapping.
#'
#' @section Aesthetics: `geom_persistence` understands the following aesthetics (required aesthetics are in bold):
#'   - **x**
#'   - **y**
#'   - point_data
#'   - alpha
#'   - color
#'   - group
#'   - size
#'
#'   `geom_barcode` understands the following aesthetics:
#'
#'   - **x**
#'   - **y**
#'   - **xend**
#'   - **yend**
#'   - point_data
#'   - alpha
#'   - color
#'   - group
#'   - linetype
#'   - size
#'
#' @section Computed variables:
#'
#'   \describe{ \item{birth}{The diameter when the feature appears}}
#'   \describe{ \item{death}{The diameter when the feature dissapears}}
#'   \describe{ \item{feature_id}{A unique id assigned to each feature}}
#'   \describe{ \item{dim}{The dimension of the feature}}
#'
#' @inheritParams ggplot2::geom_point
#' @inheritParams ggplot2::stat_identity
#' @param k Calculate features up to dimension k
#' @param point_cloud Optional, `data.frame` or `matrix` with point cloud data
#' used to calculate persitence homology. Alternatively, multiple data sets
#' can be visualized by providing mapping a list column to the `point_data`
#' aesthetic (see examples).
#'
#' @name geom_persistence
#' @rdname geom_persistence
#' @export
#'
#' @examples
#'
#' set.seed(1)
#'
#' s <- seq(0, 2*pi, length.out = 40)
#' df <- data.frame(
#'   x = cos(s) + rnorm(40, 0, .1),
#'   y = sin(s) + rnorm(40, 0, .1)
#' )
#'
#' ggplot() +
#'   geom_barcode(point_cloud = df) +
#'   scale_color_viridis_d(end = .6)
#'
#' ggplot() +
#'   geom_persistence(point_cloud = df) +
#'   geom_abline(intercept = 0, slope = 1) +
#'   coord_fixed() +
#'   scale_color_viridis_d(end = .6)
#'
#'
#' # Visualizing multiple groups together
#' # making use of point_data aesthetic mapping
#' s <- c(s, s)
#' df_mix <- data.frame(
#'   x = cos(s) + rnorm(80, 0, .1),
#'   y = sin(s) + rnorm(80, 0, .1)
#' )
#'
#' df_mix$x <- df_mix$x + rep(c(-2, 2), length.out = 80)
#' df_mix$lab <- rep(c("a", "b"), length.out = 80)
#'
#' df_nested <- dplyr::group_by(df_mix, lab)
#' df_nested <- tidyr::nest(df_nested)
#'
#' ggplot(df_nested) +
#'   geom_barcode(aes(point_data = data)) +
#'   facet_wrap(vars(lab)) +
#'   scale_color_viridis_d(end = .6)
#'
#' ggplot(df_nested) +
#'   geom_persistence(aes(point_data = data)) +
#'   geom_abline(intercept = 0, slope = 1) +
#'   coord_fixed() +
#'   facet_wrap(vars(lab)) +
#'   scale_color_viridis_d(end = .6)
#'
#' ggplot(df_nested) +
#'   geom_persistence(aes(point_data = data, shape = lab), size = 3) +
#'   geom_abline(intercept = 0, slope = 1) +
#'   coord_fixed() +
#'   scale_color_viridis_d(end = .6)
NULL

#' @rdname geom_persistence
#' @format NULL
#' @usage NULL
#' @export
StatPersistence <-  ggproto(
  "StatPersistence", Stat,

  default_aes = aes(x = after_stat(birth), y = after_stat(death), color = after_stat(dim)),

  optional_aes = c("point_data"),

  compute_group = function(data, scales, point_cloud = NULL, radius_max = NULL, diameter_max = NULL, k = 1, complex = "Rips") {

    # handle disk size
    if (is.null(radius_max) && is.null(diameter_max)) {
      diameter_max <- -1
    }
    if (! is.null(radius_max)) {
      if (! is.null(diameter_max)) {
        warning("Pass a value to only one of `radius_max` or `diameter_max`; ",
                "`diameter_max` value will be used.")
      } else {
        diameter_max <- radius_max
      }
    }

    # check if data was provided via point_cloud argument
    if (!is.null(point_cloud)){

      if (!is.data.frame(point_cloud) & !is.matrix(point_cloud)) {
        stop("Data provided to point_cloud argument must be either data.frame or matrix")
      }

      if (!is.null(data$point_data)) {
        warning(
          "Data provided via point_data aesthetic mapping and point_cloud argument.
           Defaulting to plotting point_cloud data"
        )
      }

      point_data <- point_cloud

    } else {

      if (is.null(data$point_data)) {
        stop("Must provide point data via either point_data aesthetic mapping or point_cloud argument.")
      }

      if (nrow(data) > 1) {
        # TODO - fix this message
        stop("Each group should only have one row.
             Did you forget a facet_*()?")
      }

      point_data <- data$point_data[[1]]

      if (!is.data.frame(point_data) & ! is.matrix(point_cloud)) {
        stop("Object mapped to point_data must be data.frame or matrix")
      }

    }

    res <- data_to_persistent_homology(point_data, diameter_max, k, complex)

    res <- res[order(res$dim, res$birth, res$death),]

    res$dim <- with(res, ordered(dim, levels = min(dim):max(dim)))
    res$feature_id <- 1:nrow(res)

    # returning res along with original input data
    # (recycling values from data, which only has one row!)
    rownames(data) <- NULL
    cbind(data, res)

  }
)


#' @rdname geom_persistence
#' @export
stat_persistence <- function(mapping = NULL, data = NULL,
                             geom = "persistence", position = "identity",
                             ...,
                             k = 1,
                             point_cloud = NULL,
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {

  if (is.null(data)) data <- ensure_nonempty_data

  layer(
    data = data,
    mapping = mapping,
    stat = StatPersistence,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      k = k,
      point_cloud = point_cloud,
      na.rm = na.rm,
      ...
    )
  )
}


#' @rdname geom_persistence
#' @format NULL
#' @usage NULL
#' @export
GeomPersistence <-  ggproto(
  "GeomPersistence", GeomPoint,
  required_aes = c("x", "y"),
  optional_aes = c("point_data")
)


#' @rdname geom_persistence
#' @export
geom_persistence <- function(mapping = NULL, data = NULL,
                             stat = "persistence", position = "identity",
                             ...,
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {

  if (is.null(data)) data <- ensure_nonempty_data

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomPersistence,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}




# Barcode charts ----------------------------------------------------------
# inherits computations from StatPersistence

#' @rdname geom_persistence
#' @format NULL
#' @usage NULL
#' @export
StatBarcode <-  ggproto(
  "StatBarcode", StatPersistence,

  default_aes = aes(x = after_stat(birth), xend = after_stat(death),
                    y = after_stat(feature_id), yend = after_stat(feature_id),
                    color = after_stat(dim))

)


#' @rdname geom_persistence
#' @export
stat_barcode <- function(mapping = NULL, data = NULL,
                         geom = "barcode", position = "identity",
                         ...,
                         k = 1,
                         point_cloud = NULL,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {

  if (is.null(data)) data <- ensure_nonempty_data

  layer(
    data = data,
    mapping = mapping,
    stat = StatBarcode,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      k = k,
      point_cloud = point_cloud,
      na.rm = na.rm,
      ...
    )
  )
}


#' @rdname geom_persistence
#' @format NULL
#' @usage NULL
#' @export
GeomBarcode <-  ggproto(
  "GeomBarcode", GeomSegment,
  required_aes = c("x", "y", "xend", "yend"),
  optional_aes = c("point_data")
)


#' @rdname geom_persistence
#' @export
geom_barcode <- function(mapping = NULL, data = NULL,
                         stat = "barcode", position = "identity",
                         ...,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {

  if (is.null(data)) data <- ensure_nonempty_data

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomBarcode,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

