#' Title
#'
#' Compute persistence homologies and plot barcode plots.
#' Accepts either a singular point cloud of arbitrary dimension via the `point_cloud` argument
#' or several, via the `point_data` aesthetic mapping.
#'
#' @section Aesthetics: geom_hdr understands the following aesthetics (required
#'   aesthetics are in bold):
#'
#'   - point_data
#'   - alpha
#'   - color
#'   - group
#'   - linetype
#'   - size
#'   - x
#'   - y
#'   - xend
#'   - yend
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
#' @param k
#' @param point_cloud
#'
#' @name geom_barcode
#' @rdname geom_barcode
#' @export
#'
#' @examples
#'
#' # generate a noisy circle
#' n <- 36; sd <- .2
#' set.seed(0)
#' t <- stats::runif(n = n, min = 0, max = 2*pi)
#' df <- data.frame(
#'   x = cos(t) + stats::rnorm(n = n, mean = 0, sd = sd),
#'   y = sin(t) + stats::rnorm(n = n, mean = 0, sd = sd)
#' )
#'
#' ggplot() +
#'   geom_barcode(point_cloud = df) +
#'   scale_color_viridis_d(end = .6)
#'
#'
#' # Visualizing multiple groups together
#' # making use of point_data aesthetic mapping
#' n <- 72; sd <- .2
#' set.seed(0)
#' t <- stats::runif(n = n, min = 0, max = 2*pi)
#' df_mix <- data.frame(
#'   x = cos(t) + stats::rnorm(n = n, mean = 0, sd = sd),
#'   y = sin(t) + stats::rnorm(n = n, mean = 0, sd = sd)
#' )
#'
#' df_mix$x <- df_mix$x + rep(c(-2, 2), length.out = n)
#' df_mix$lab <- rep(c("a", "b"), length.out = n)
#'
#' df_nested <- dplyr::group_by(df_mix, lab)
#' df_nested <- tidyr::nest(df_nested)
#'
#' ggplot() +
#'   geom_barcode(data = df_nested, aes(point_data = data)) +
#'   facet_wrap(vars(lab)) +
#'   scale_color_viridis_d(end = .6)
NULL

#' @rdname geom_barcode
#' @export
StatBarcode <-  ggproto(
  "StatBarcode", Stat,

  default_aes = aes(x = after_stat(birth), xend = after_stat(death),
                    y = after_stat(feature_id), yend = after_stat(feature_id),
                    color = after_stat(dim)),

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


#' @rdname geom_barcode
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


#' @rdname geom_barcode
#' @export
GeomBarcode <-  ggproto(
  "GeomBarcode", GeomSegment,
  optional_aes = c("point_data")
)


#' @rdname geom_barcode
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
