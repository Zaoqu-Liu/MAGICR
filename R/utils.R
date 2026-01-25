#' @title MAGICR Utility Functions
#' @description Internal utility functions for the MAGICR package.
#' @name utils
#' @keywords internal
NULL

#' Check if a value is positive
#' @param x Value to check
#' @param name Parameter name for error message
#' @keywords internal
check_positive <- function(x, name = deparse(substitute(x))) {
  if (!is.null(x) && (length(x) != 1 || !is.numeric(x) || x <= 0)) {
    stop(sprintf("'%s' must be a positive number, got: %s", name, x),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check if a value is a positive integer
#' @param x Value to check
#' @param name Parameter name for error message
#' @keywords internal
check_positive_int <- function(x, name = deparse(substitute(x))) {
  if (!is.null(x) && (length(x) != 1 || !is.numeric(x) || x <= 0 || x != floor(x))) {
    stop(sprintf("'%s' must be a positive integer, got: %s", name, x),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check if a value is in a set of choices
#' @param x Value to check
#' @param choices Valid choices
#' @param name Parameter name for error message
#' @keywords internal
check_in <- function(x, choices, name = deparse(substitute(x))) {
  if (!x %in% choices) {
    stop(sprintf("'%s' must be one of: %s. Got: %s",
                 name, paste(choices, collapse = ", "), x),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Safe conversion to sparse matrix
#' @param x Matrix or data.frame
#' @return dgCMatrix
#' @keywords internal
to_sparse <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  if (inherits(x, "Matrix")) {
    return(as(x, "dgCMatrix"))
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  return(Matrix::Matrix(x, sparse = TRUE))
}

#' Safe conversion to dense matrix
#' @param x Matrix or data.frame
#' @return matrix
#' @keywords internal
to_dense <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }
  if (inherits(x, "Matrix")) {
    return(as.matrix(x))
  }
  if (is.data.frame(x)) {
    return(as.matrix(x))
  }
  return(as.matrix(x))
}

#' Check if matrix has empty columns (unexpressed genes)
#' @param x Matrix
#' @return logical
#' @keywords internal
has_empty_columns <- function(x) {
  col_sums <- Matrix::colSums(x)
  any(col_sums == 0)
}

#' Check if matrix has empty rows (cells with no counts)
#' @param x Matrix
#' @return logical
#' @keywords internal
has_empty_rows <- function(x) {
  row_sums <- Matrix::rowSums(x)
  any(row_sums == 0)
}

#' Get row names safely
#' @param x Matrix or data.frame
#' @return character vector or NULL
#' @keywords internal
safe_rownames <- function(x) {
  rn <- rownames(x)
  if (is.null(rn) || length(rn) == 0) {
    return(paste0("Cell_", seq_len(nrow(x))))
  }
  return(rn)
}

#' Get column names safely
#' @param x Matrix or data.frame
#' @return character vector or NULL
#' @keywords internal
safe_colnames <- function(x) {
  cn <- colnames(x)
  if (is.null(cn) || length(cn) == 0) {
    return(paste0("Gene_", seq_len(ncol(x))))
  }
  return(cn)
}

#' Verbose message printing
#' @param msg Message to print
#' @param verbose Verbosity level
#' @param level Required level to print
#' @keywords internal
vmessage <- function(msg, verbose = 1, level = 1) {
  if (is.logical(verbose)) {
    verbose <- as.integer(verbose)
  }
  if (verbose >= level) {
    message(msg)
  }
  invisible(NULL)
}

#' Timer utility for performance logging
#' @param expr Expression to time
#' @param name Name for the operation
#' @param verbose Verbosity level
#' @return Result of expression
#' @keywords internal
timed <- function(expr, name = "Operation", verbose = 1) {
  if (verbose >= 2) {
    start_time <- Sys.time()
    result <- force(expr)
    elapsed <- difftime(Sys.time(), start_time, units = "secs")
    message(sprintf("  [%s] completed in %.2f seconds", name, as.numeric(elapsed)))
    return(result)
  } else {
    return(force(expr))
  }
}

#' Check if Rcpp functions are available
#' @return logical
#' @keywords internal
rcpp_available <- function() {
  tryCatch({
    .Call("_MAGICR_rcpp_available", PACKAGE = "MAGICR")
    TRUE
  }, error = function(e) FALSE)
}

#' Null-coalescing operator
#' @param x First value
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Check if parallel backend is configured
#' @return logical
#' @keywords internal
is_parallel <- function() {
  !inherits(future::plan(), "sequential")
}

#' Get number of workers
#' @return integer
#' @keywords internal
n_workers <- function() {
  future::nbrOfWorkers()
}

#' Procrustes disparity calculation
#' @description Calculate the Procrustes disparity between two matrices.
#'   Used for automatic t selection in MAGIC.
#' @param X First matrix (n x p)
#' @param Y Second matrix (n x p)
#' @return Procrustes disparity value (0 to 1)
#' @keywords internal
procrustes_disparity <- function(X, Y) {
  # Center both matrices
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  
  # Compute Frobenius norms
  X_norm <- sqrt(sum(X^2))
  Y_norm <- sqrt(sum(Y^2))
  
  # Normalize
  X <- X / X_norm
  Y <- Y / Y_norm
  

  # SVD of X'Y
  svd_result <- svd(crossprod(X, Y))
  
  # Optimal rotation
  R <- tcrossprod(svd_result$v, svd_result$u)
  
  # Disparity: 1 - trace(R * svd_result$d) but simplified
  # disparity = sum((X - Y %*% R)^2) after alignment
  Y_aligned <- Y %*% R
  disparity <- sum((X - Y_aligned)^2)
  
  return(disparity)
}

#' Select columns by index or name
#' @param x Matrix or data.frame
#' @param idx Integer indices or character names
#' @return Subsetted matrix
#' @keywords internal
select_cols <- function(x, idx) {
  if (is.null(idx)) {
    return(x)
  }
  if (is.character(idx)) {
    if (!all(idx %in% colnames(x))) {
      missing <- idx[!idx %in% colnames(x)]
      warning(sprintf("Genes not found: %s", paste(missing, collapse = ", ")))
      idx <- idx[idx %in% colnames(x)]
    }
  }
  x[, idx, drop = FALSE]
}

#' Platform-independent memory check
#' @param bytes Number of bytes needed
#' @return logical indicating if memory is likely available
#' @keywords internal
check_memory <- function(bytes) {
  # Conservative estimate - just warn for very large allocations
  gb <- bytes / (1024^3)
  if (gb > 8) {
    warning(sprintf("This operation may require approximately %.1f GB of memory.", gb))
  }
  invisible(TRUE)
}

#' Format time duration
#' @param seconds Number of seconds
#' @return Formatted string
#' @keywords internal
format_time <- function(seconds) {
  if (seconds < 60) {
    sprintf("%.1f seconds", seconds)
  } else if (seconds < 3600) {
    sprintf("%.1f minutes", seconds / 60)
  } else {
    sprintf("%.1f hours", seconds / 3600)
  }
}
