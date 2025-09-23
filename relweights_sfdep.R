#' Compute Relational Weights using sfdep
#'
#' This implementation constructs relational weights and returns structures that
#' integrate with the `sfdep` workflow. A custom GAL writer is provided so users can
#' persist the neighbour graph without relying on `spdep`.
#'
#' @param rel_layer sf object; the layer inheriting neighbours.
#' @param inh_layer sf object; the layer inducing inherited connections.
#' @param style Character; currently supports "W" (row-standardised) and "B" (binary).
#' @param binary Logical; when TRUE (default) cell values greater than zero collapse to 1 prior to
#'   applying the chosen style.
#' @param zero.policy Logical; when FALSE the presence of zero-neighbour observations triggers an
#'   error when style = "W" and a row standardisation would divide by zero.
#' @param ids Optional character vector of identifiers to label the rows when exporting.
#' @param output_gal Optional path; when provided a GAL file is written using the computed
#'   neighbour structure.
#'
#' @return A list containing the dense relational matrix, the raw weights matrix, the row-standardised
#'   (or binary) weights matrix used to build sfdep objects, the neighbour list with class `nb`, a list of
#'   weights compatible with sfdep, a recreated listw object, and a data.frame with id/nb/weight columns
#'   for tidy workflows.
#' @export
#'
#' @examples
#' # result <- tab2relweights_sfdep(rel_sf, inh_sf, output_gal = "RelWeights.gal")
#' # result$nb  # sfdep-compatible neighbour list

tab2relweights_sfdep <- function(rel_layer,
                                 inh_layer,
                                 style = c("W", "B"),
                                 binary = TRUE,
                                 zero.policy = TRUE,
                                 ids = NULL,
                                 output_gal = NULL) {
  style <- match.arg(style)

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.")
  }
  if (!requireNamespace("sfdep", quietly = TRUE)) {
    stop("Package 'sfdep' is required.")
  }
  if (!inherits(rel_layer, "sf")) {
    stop("`rel_layer` must be an sf object.")
  }
  if (!inherits(inh_layer, "sf")) {
    stop("`inh_layer` must be an sf object.")
  }

  n_rel <- nrow(rel_layer)
  n_inh <- nrow(inh_layer)

  if (!is.null(ids)) {
    if (length(ids) != n_rel) {
      stop("`ids` must have length equal to nrow(rel_layer).")
    }
    rel_ids <- as.character(ids)
  } else {
    rel_ids <- as.character(seq_len(n_rel))
  }

  rel_layer$id_x <- seq_len(n_rel)
  inh_layer$id_y <- seq_len(n_inh)

  hits <- sf::st_intersects(rel_layer, inh_layer)

  rel_matrix <- matrix(0, nrow = n_rel, ncol = n_rel)
  base_weights <- rel_matrix
  nb_list <- replicate(n_rel, integer(0), simplify = FALSE)

  if (!all(lengths(hits) == 0)) {
    incidence <- build_incidence_matrix(n_rel, n_inh, hits)

    rel_matrix <- incidence %*% t(incidence)
    diag(rel_matrix) <- 0

    nb_list <- lapply(seq_len(n_rel), function(i) which(rel_matrix[i, ] > 0))

    if (binary) {
      base_weights <- ifelse(rel_matrix > 0, 1, 0)
    } else {
      base_weights <- rel_matrix
    }
  }

  styled_weights <- apply_style(base_weights, style, zero.policy)
  nb_list <- ensure_nb(nb_list, rel_ids)
  weights_list <- matrix_to_weight_list(styled_weights, nb_list, style)

  listw <- sfdep::recreate_listw(nb_list, weights_list)
  tidy_df <- data.frame(
    id = rel_ids,
    nb = I(nb_list),
    wt = I(weights_list),
    stringsAsFactors = FALSE
  )

  if (!is.null(output_gal)) {
    write_gal(nb_list, rel_ids, output_gal)
  }

  list(rel_matrix = rel_matrix,
       weights = base_weights,
       weights_styled = styled_weights,
       nb = nb_list,
       weights_list = weights_list,
       listw = listw,
       sfdep = tidy_df)
}

# Helper: ensure neighbour list carries expected metadata
ensure_nb <- function(nb_list, ids) {
  nb_list <- lapply(nb_list, function(x) as.integer(x))
  class(nb_list) <- c("nb", "list")
  attr(nb_list, "region.id") <- ids
  attr(nb_list, "sym") <- TRUE
  attr(nb_list, "type") <- "relational"
  nb_list
}

# Helper: convert weights matrix to sfdep-compatible list representation
matrix_to_weight_list <- function(mat, nb_list, style) {
  n <- nrow(mat)
  wt <- vector("list", length = n)
  for (i in seq_len(n)) {
    neighbours <- nb_list[[i]]
    if (length(neighbours) == 0) {
      wt[[i]] <- numeric(0)
    } else {
      wt[[i]] <- mat[i, neighbours]
    }
  }
  attr(wt, style) <- TRUE
  wt
}

# Helper: apply weighting style with basic support for W and B
apply_style <- function(mat, style, zero.policy) {
  if (style == "B") {
    mat[mat > 0] <- 1
    return(mat)
  }
  if (style != "W") {
    stop("Currently only styles 'W' and 'B' are supported in the sfdep implementation.")
  }
  row_sums <- rowSums(mat)
  zero_rows <- row_sums == 0
  if (any(zero_rows) && !zero.policy) {
    stop("Zero neighbours encountered and zero.policy = FALSE prevents row standardisation.")
  }
  row_sums[zero_rows] <- 1
  mat / row_sums
}

# Helper: minimal GAL writer working with nb + ids
write_gal <- function(nb_list, ids, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  for (i in seq_along(nb_list)) {
    neighbours <- nb_list[[i]]
    line_one <- sprintf("%s %d", ids[i], length(neighbours))
    writeLines(line_one, con)
    if (length(neighbours) == 0) {
      writeLines("", con)
    } else {
      neighbour_ids <- ids[neighbours]
      writeLines(paste(neighbour_ids, collapse = " "), con)
    }
  }
  invisible(TRUE)
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
