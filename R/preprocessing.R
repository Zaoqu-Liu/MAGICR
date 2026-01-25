#' @title Data Preprocessing Functions
#' @description Functions for preprocessing single-cell data before MAGIC.
#' @name preprocessing
NULL

#' Library Size Normalization
#'
#' @description Performs library size (L1) normalization on input data such that
#'   the sum of expression values for each cell sums to 1, then returns normalized
#'   matrix scaled by median library size.
#'
#' @param data A matrix (cells x genes) of expression values. Can be a dense
#'   matrix, sparse Matrix, or data.frame.
#' @param verbose Logical or integer. If TRUE or > 0, print progress messages.
#'
#' @return A matrix of the same type and dimensions as input, with normalized
#'   expression values.
#'
#' @details Library size normalization is a standard preprocessing step for
#'   single-cell RNA-seq data. It normalizes each cell by its total UMI count
#'   and then scales by the median total count across all cells. This makes
#'   cells comparable by removing technical variation in sequencing depth.
#'
#'   The formula is: \deqn{X_{norm} = \frac{X}{rowSums(X)} \times median(rowSums(X))}
#'
#' @export
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(42)
#' data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#'
#' # Normalize
#' data_norm <- library_size_normalize(data)
#'
#' # Check that median library size is preserved
#' all.equal(median(rowSums(data)), median(rowSums(data_norm)))
#' }
library_size_normalize <- function(data, verbose = FALSE) {
  vmessage("Performing library size normalization...", verbose)
  
  # Handle sparse matrices efficiently
  is_sparse <- inherits(data, "sparseMatrix")
  
  # Calculate library sizes (total counts per cell)
  library_size <- Matrix::rowSums(data)
  
  # Check for zero library sizes
  if (any(library_size == 0)) {
    warning("Some cells have zero total counts. These will result in NaN values.")
  }
  
  # Calculate median library size
  median_size <- stats::median(library_size[library_size > 0])
  
  # Normalize: divide by library size, multiply by median
  # Vectorized operation that preserves sparsity
  if (is_sparse) {
    # For sparse matrices, use efficient diagonal scaling
    scale_factors <- median_size / library_size
    scale_factors[!is.finite(scale_factors)] <- 0
    data_norm <- Matrix::Diagonal(x = scale_factors) %*% data
  } else {
    data_norm <- data * (median_size / library_size)
    data_norm[!is.finite(data_norm)] <- 0
  }
  
  # Preserve row/column names
  rownames(data_norm) <- rownames(data)
  colnames(data_norm) <- colnames(data)
  
  vmessage(sprintf("  Normalized %d cells, median library size: %.1f",
                   nrow(data), median_size), verbose)
  
  return(data_norm)
}

#' Square Root Normalization
#'
#' @description Applies square root transformation to expression data.
#'   This is a variance-stabilizing transformation commonly used for count data.
#'
#' @param data A matrix (cells x genes) of expression values.
#' @param verbose Logical or integer. If TRUE or > 0, print progress messages.
#'
#' @return A matrix of the same dimensions with square root transformed values.
#'
#' @details Square root transformation is useful for stabilizing variance in
#'   count data where the variance tends to increase with the mean. It is
#'   less aggressive than log transformation and handles zeros naturally.
#'
#' @export
#' @examples
#' \dontrun{
#' data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#' data_sqrt <- sqrt_normalize(data)
#' }
sqrt_normalize <- function(data, verbose = FALSE) {
  vmessage("Applying square root transformation...", verbose)
  
  # Check for negative values
  if (any(data < 0)) {
    warning("Negative values detected. These will result in NaN values.")
  }
  
  # Apply transformation (preserves sparsity structure for sparse matrices)
  data_sqrt <- sqrt(data)
  
  # Preserve row/column names
  rownames(data_sqrt) <- rownames(data)
  colnames(data_sqrt) <- colnames(data)
  
  return(data_sqrt)
}

#' Log Transformation
#'
#' @description Applies log transformation with a pseudo-count to expression data.
#'
#' @param data A matrix (cells x genes) of expression values.
#' @param pseudo_count Numeric. Value added before log transformation to avoid
#'   log(0). Default is 1.
#' @param base Numeric. Base of the logarithm. Default is natural log (exp(1)).
#' @param verbose Logical or integer. If TRUE or > 0, print progress messages.
#'
#' @return A matrix of the same dimensions with log transformed values.
#'
#' @details Log transformation is commonly used to normalize the distribution
#'   of gene expression values. The pseudo-count ensures that zeros don't
#'   produce -Inf values.
#'
#'   Formula: \deqn{X_{log} = \log_{base}(X + pseudo\_count)}
#'
#' @export
#' @examples
#' \dontrun{
#' data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#' data_log <- log_transform(data, pseudo_count = 1)
#' }
log_transform <- function(data, pseudo_count = 1, base = exp(1), verbose = FALSE) {
  vmessage(sprintf("Applying log transformation (pseudo_count = %g, base = %g)...",
                   pseudo_count, base), verbose)
  
  check_positive(pseudo_count, "pseudo_count")
  check_positive(base, "base")
  
  # Check for negative values
  if (any(data < 0)) {
    warning("Negative values detected. Adding pseudo_count may not be sufficient.")
  }
  
  # Handle sparse matrices - log transform destroys sparsity
  if (inherits(data, "sparseMatrix")) {
    vmessage("  Note: Log transformation converts sparse matrix to dense.", verbose)
    data <- as.matrix(data)
  }
  
  # Apply transformation
  data_log <- log(data + pseudo_count, base = base)
  
  # Preserve row/column names
  rownames(data_log) <- rownames(data)
  colnames(data_log) <- colnames(data)
  
  return(data_log)
}

#' Remove Low-Count Genes
#'
#' @description Filters out genes expressed in fewer than a minimum number of cells.
#'
#' @param data A matrix (cells x genes) of expression values.
#' @param min_cells Integer. Minimum number of cells in which a gene must be
#'   expressed (count > 0). Default is 10.
#' @param verbose Logical or integer. If TRUE or > 0, print progress messages.
#'
#' @return A matrix with low-count genes removed.
#'
#' @keywords internal
filter_genes <- function(data, min_cells = 10, verbose = FALSE) {
  # Count cells with non-zero expression for each gene
  if (inherits(data, "sparseMatrix")) {
    n_cells <- Matrix::colSums(data > 0)
  } else {
    n_cells <- colSums(data > 0)
  }
  
  keep <- n_cells >= min_cells
  n_removed <- sum(!keep)
  
  if (n_removed > 0) {
    vmessage(sprintf("  Removed %d genes expressed in fewer than %d cells",
                     n_removed, min_cells), verbose)
    data <- data[, keep, drop = FALSE]
  }
  
  return(data)
}

#' Remove Low-Count Cells
#'
#' @description Filters out cells with fewer than a minimum number of detected genes.
#'
#' @param data A matrix (cells x genes) of expression values.
#' @param min_genes Integer. Minimum number of genes that must be detected
#'   (count > 0) in a cell. Default is 200.
#' @param verbose Logical or integer. If TRUE or > 0, print progress messages.
#'
#' @return A matrix with low-count cells removed.
#'
#' @keywords internal
filter_cells <- function(data, min_genes = 200, verbose = FALSE) {
  # Count genes with non-zero expression for each cell
  if (inherits(data, "sparseMatrix")) {
    n_genes <- Matrix::rowSums(data > 0)
  } else {
    n_genes <- rowSums(data > 0)
  }
  
  keep <- n_genes >= min_genes
  n_removed <- sum(!keep)
  
  if (n_removed > 0) {
    vmessage(sprintf("  Removed %d cells with fewer than %d genes detected",
                     n_removed, min_genes), verbose)
    data <- data[keep, , drop = FALSE]
  }
  
  return(data)
}

#' Validate Input Data
#'
#' @description Validates and prepares input data for MAGIC.
#'
#' @param data Input data matrix
#' @param verbose Verbosity level
#'
#' @return Validated matrix
#' @keywords internal
validate_data <- function(data, verbose = FALSE) {
  # Convert to matrix if needed
  if (is.data.frame(data)) {
    # Check for non-numeric columns
    numeric_cols <- vapply(data, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      vmessage("  Converting non-numeric columns to numeric", verbose)
      data <- data[, numeric_cols, drop = FALSE]
    }
    data <- as.matrix(data)
  }
  
  # Check dimensions
  if (nrow(data) < 10) {
    warning("Input data has fewer than 10 cells. MAGIC may not work well.")
  }
  if (ncol(data) < 10) {
    warning("Input data has fewer than 10 genes. MAGIC may not work well.")
  }
  
  # Check for NA/NaN/Inf
  if (any(!is.finite(data))) {
    n_bad <- sum(!is.finite(data))
    warning(sprintf("Input contains %d non-finite values (NA/NaN/Inf). Setting to 0.",
                    n_bad))
    data[!is.finite(data)] <- 0
  }
  
  # Check for negative values
  if (any(data < 0)) {
    warning("Input contains negative values. This may indicate already normalized data.")
  }
  
  # Check for empty columns
  if (has_empty_columns(data)) {
    warning("Input contains unexpressed genes (all-zero columns). ",
            "Consider removing them before running MAGIC.")
  }
  
  # Check for empty rows
  if (has_empty_rows(data)) {
    warning("Input contains cells with no counts (all-zero rows). ",
            "Consider removing them before running MAGIC.")
  }
  
  return(data)
}
