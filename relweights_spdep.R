#' Compute Relational Weights using spdep
#'
#' This implementation uses `spdep` to construct listw objects and optionally
#' export old-style GAL files.
#'
#' @param rel_layer sf object representing the layer that inherits neighbours.
#' @param inh_layer sf object representing the layer that induces connections.
#' @param style Character, weighting style forwarded to `spdep::mat2listw`.
#' @param binary Logical; when TRUE (default) weights greater than zero become 1.
#' @param zero.policy Logical; passed to `spdep::mat2listw` for zero-neighbour rows.
#' @param output_gal Optional path; when supplied a GAL file is written.
#'
#' @return A list with the relational matrix, the weights matrix, the `spdep`
#'   listw object, and the neighbour list.
#' @export
#'
#' @examples
#' # result <- tab2relweights_spdep(rel_sf, inh_sf, output_gal = "RelWeights.gal")
#' # result$listw

tab2relweights_spdep <- function(rel_layer,
                                 inh_layer,
                                 style = "W",
                                 binary = TRUE,
                                 zero.policy = TRUE,
                                 output_gal = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.")
  }
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required.")
  }
  if (!inherits(rel_layer, "sf")) {
    stop("`rel_layer` must be an sf object.")
  }
  if (!inherits(inh_layer, "sf")) {
    stop("`inh_layer` must be an sf object.")
  }

  n_rel <- nrow(rel_layer)
  n_inh <- nrow(inh_layer)

  rel_layer$id_x <- seq_len(n_rel)
  inh_layer$id_y <- seq_len(n_inh)

  hits <- sf::st_intersects(rel_layer, inh_layer)

  if (all(lengths(hits) == 0)) {
    zero_mat <- matrix(0, nrow = n_rel, ncol = n_rel)
    listw <- spdep::mat2listw(zero_mat, style = style, zero.policy = zero.policy)
    if (!is.null(output_gal)) {
      spdep::write.nb.gal(listw$neighbours, file = output_gal, oldstyle = TRUE)
    }
    return(list(rel_matrix = zero_mat,
                weights = zero_mat,
                listw = listw,
                nb = listw$neighbours))
  }

  incidence <- build_incidence_matrix(n_rel, n_inh, hits)

  rel_matrix <- incidence %*% t(incidence)
  diag(rel_matrix) <- 0

  if (binary) {
    weights_matrix <- ifelse(rel_matrix > 0, 1, 0)
  } else {
    weights_matrix <- rel_matrix
  }

  listw <- spdep::mat2listw(weights_matrix, style = style, zero.policy = zero.policy)

  if (!is.null(output_gal)) {
    spdep::write.nb.gal(listw$neighbours, file = output_gal, oldstyle = TRUE)
  }

  list(rel_matrix = rel_matrix,
       weights = weights_matrix,
       listw = listw,
       nb = listw$neighbours)
}

build_incidence_matrix <- function(n_rel, n_inh, hits) {
  incidence <- matrix(0, nrow = n_rel, ncol = n_inh)
  idx_rel <- rep.int(seq_along(hits), lengths(hits))
  idx_inh <- unlist(hits, use.names = FALSE)
  if (length(idx_rel) > 0) {
    indices <- cbind(idx_rel, idx_inh)
    incidence[indices] <- incidence[indices] + 1
  }
  incidence
}
