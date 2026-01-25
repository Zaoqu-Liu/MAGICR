#' @title Data Preprocessing
#' @description Preprocessing functions for single-cell data.
#' @name preprocessing
NULL

#' Library Size Normalization
#'
#' @description Normalizes data by library size (total counts per cell).
#'
#' @param data Matrix (cells x genes).
#' @param verbose Logical. Print progress. Default is FALSE.
#'
#' @return Normalized matrix.
#'
#' @details Formula: \eqn{X_{norm} = X / rowSums(X) * median(rowSums(X))}
#'
#' @export
library_size_normalize <- function(data, verbose = FALSE) {
  vmessage("Library size normalization...", verbose)
  
  is_sparse <- inherits(data, "sparseMatrix")
  
  library_size <- Matrix::rowSums(data)
  
  if (any(library_size == 0)) {
    warning("Some cells have zero counts.")
  }
  
  median_size <- stats::median(library_size[library_size > 0])
  
  if (is_sparse) {
    scale_factors <- median_size / library_size
    scale_factors[!is.finite(scale_factors)] <- 0
    data_norm <- Matrix::Diagonal(x = scale_factors) %*% data
  } else {
    data_norm <- data * (median_size / library_size)
    data_norm[!is.finite(data_norm)] <- 0
  }
  
  rownames(data_norm) <- rownames(data)
  colnames(data_norm) <- colnames(data)
  
  vmessage(sprintf("  Normalized %d cells, median size: %.1f",
                   nrow(data), median_size), verbose)
  
  return(data_norm)
}

#' Square Root Transformation
#'
#' @description Applies square root transformation.
#'
#' @param data Matrix.
#' @param verbose Logical. Default is FALSE.
#'
#' @return Transformed matrix.
#'
#' @export
sqrt_normalize <- function(data, verbose = FALSE) {
  vmessage("Square root transformation...", verbose)
  
  if (any(data < 0)) {
    warning("Negative values detected.")
  }
  
  data_sqrt <- sqrt(data)
  rownames(data_sqrt) <- rownames(data)
  colnames(data_sqrt) <- colnames(data)
  
  return(data_sqrt)
}

#' Log Transformation
#'
#' @description Applies log transformation with pseudo-count.
#'
#' @param data Matrix.
#' @param pseudo_count Numeric. Default is 1.
#' @param base Numeric. Log base. Default is exp(1).
#' @param verbose Logical. Default is FALSE.
#'
#' @return Transformed matrix.
#'
#' @details Formula: \eqn{X_{log} = log_{base}(X + pseudo\_count)}
#'
#' @export
log_transform <- function(data, pseudo_count = 1, base = exp(1), verbose = FALSE) {
  vmessage(sprintf("Log transformation (pseudo=%g)...", pseudo_count), verbose)
  
  check_positive(pseudo_count, "pseudo_count")
  check_positive(base, "base")
  
  if (any(data < 0)) {
    warning("Negative values detected.")
  }
  
  if (inherits(data, "sparseMatrix")) {
    vmessage("  Converting sparse to dense for log transform.", verbose)
    data <- as.matrix(data)
  }
  
  data_log <- log(data + pseudo_count, base = base)
  rownames(data_log) <- rownames(data)
  colnames(data_log) <- colnames(data)
  
  return(data_log)
}

#' Filter genes by expression
#' @param data Matrix
#' @param min_cells Minimum cells. Default is 10.
#' @param verbose Logical. Default is FALSE.
#' @return Filtered matrix
#' @keywords internal
filter_genes <- function(data, min_cells = 10, verbose = FALSE) {
  if (inherits(data, "sparseMatrix")) {
    n_cells <- Matrix::colSums(data > 0)
  } else {
    n_cells <- colSums(data > 0)
  }
  
  keep <- n_cells >= min_cells
  n_removed <- sum(!keep)
  
  if (n_removed > 0) {
    vmessage(sprintf("  Removed %d genes (< %d cells)", n_removed, min_cells), verbose)
    data <- data[, keep, drop = FALSE]
  }
  
  return(data)
}

#' Filter cells by gene count
#' @param data Matrix
#' @param min_genes Minimum genes. Default is 200.
#' @param verbose Logical. Default is FALSE.
#' @return Filtered matrix
#' @keywords internal
filter_cells <- function(data, min_genes = 200, verbose = FALSE) {
  if (inherits(data, "sparseMatrix")) {
    n_genes <- Matrix::rowSums(data > 0)
  } else {
    n_genes <- rowSums(data > 0)
  }
  
  keep <- n_genes >= min_genes
  n_removed <- sum(!keep)
  
  if (n_removed > 0) {
    vmessage(sprintf("  Removed %d cells (< %d genes)", n_removed, min_genes), verbose)
    data <- data[keep, , drop = FALSE]
  }
  
  return(data)
}

#' Validate input data
#' @param data Input data
#' @param verbose Verbosity level
#' @return Validated matrix
#' @keywords internal
validate_data <- function(data, verbose = FALSE) {
  if (is.data.frame(data)) {
    numeric_cols <- vapply(data, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      data <- data[, numeric_cols, drop = FALSE]
    }
    data <- as.matrix(data)
  }
  
  if (nrow(data) < 10) {
    warning("Fewer than 10 cells.")
  }
  if (ncol(data) < 10) {
    warning("Fewer than 10 genes.")
  }
  
  if (any(!is.finite(data))) {
    n_bad <- sum(!is.finite(data))
    warning(sprintf("Setting %d non-finite values to 0.", n_bad))
    data[!is.finite(data)] <- 0
  }
  
  if (any(data < 0)) {
    warning("Negative values detected.")
  }
  
  if (has_empty_columns(data)) {
    warning("Some genes have zero expression.")
  }
  
  if (has_empty_rows(data)) {
    warning("Some cells have zero counts.")
  }
  
  return(data)
}
