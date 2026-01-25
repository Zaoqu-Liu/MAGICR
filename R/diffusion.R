#' @title Diffusion Operations
#' @description Functions for diffusion-based data imputation.
#' @name diffusion
NULL

#' Perform MAGIC Imputation
#'
#' @description Applies the diffusion operator to impute and denoise data.
#'
#' @param data Matrix (cells x genes) to impute.
#' @param diff_op Diffusion operator (Markov matrix).
#' @param t Integer or "auto". Number of diffusion steps.
#' @param t_max Integer. Maximum t for automatic selection. Default is 20.
#' @param threshold Numeric. Convergence threshold for automatic t. Default is 0.001.
#' @param solver Character. "exact" or "approximate". Default is "exact".
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A matrix of imputed values.
#'
#' @details MAGIC imputation applies the diffusion operator repeatedly:
#'   \deqn{X_{imputed} = P^t \times X}
#'   
#'   When t = "auto", the algorithm stops when the Procrustes disparity
#'   between consecutive iterations falls below the threshold.
#'
#' @export
magic_impute <- function(data,
                         diff_op,
                         t = 3,
                         t_max = 20,
                         threshold = 0.001,
                         solver = c("exact", "approximate"),
                         verbose = FALSE) {
  
  solver <- match.arg(solver)
  
  # Validate inputs
  if (!inherits(diff_op, "Matrix") && !is.matrix(diff_op)) {
    stop("diff_op must be a matrix", call. = FALSE)
  }
  
  if (nrow(diff_op) != ncol(diff_op)) {
    stop("diff_op must be square", call. = FALSE)
  }
  
  if (nrow(diff_op) != nrow(data)) {
    stop("Number of rows in data must match diff_op dimensions", call. = FALSE)
  }
  
  # Handle t parameter
  auto_t <- identical(t, "auto")
  if (!auto_t) {
    check_positive_int(t, "t")
  }
  check_positive_int(t_max, "t_max")
  check_positive(threshold, "threshold")
  
  # Convert to dense
  data <- to_dense(data)
  
  if (auto_t) {
    result <- impute_with_auto_t(data, diff_op, t_max, threshold, verbose)
  } else {
    result <- impute_with_t(data, diff_op, t, verbose)
  }
  
  return(result)
}

#' Imputation with specified t
#' @param data Data matrix
#' @param diff_op Diffusion operator
#' @param t Number of steps
#' @param verbose Verbosity level
#' @return Imputed matrix
#' @keywords internal
impute_with_t <- function(data, diff_op, t, verbose = FALSE) {
  n_cells <- nrow(data)
  n_genes <- ncol(data)
  
  vmessage(sprintf("Imputing with t=%d...", t), verbose)
  
  # Choose strategy based on matrix dimensions
  if (n_cells < n_genes && t > 1) {
    # Power the operator first, then multiply
    P_t <- matrix_power(diff_op, t)
    result <- as.matrix(P_t %*% data)
  } else {
    # Iterative multiplication
    result <- data
    for (i in seq_len(t)) {
      result <- as.matrix(diff_op %*% result)
    }
  }
  
  # Preserve names
  rownames(result) <- rownames(data)
  colnames(result) <- colnames(data)
  
  vmessage("  Imputation complete", verbose)
  
  return(result)
}

#' Imputation with automatic t selection
#' @param data Data matrix
#' @param diff_op Diffusion operator
#' @param t_max Maximum t
#' @param threshold Convergence threshold
#' @param verbose Verbosity level
#' @return Imputed matrix
#' @keywords internal
impute_with_auto_t <- function(data, diff_op, t_max, threshold, verbose = FALSE) {
  n_genes <- ncol(data)
  
  vmessage("Selecting optimal t automatically...", verbose)
  
  # Subsample genes for faster error calculation
  max_genes <- 500
  if (n_genes > max_genes) {
    subsample_idx <- sample(n_genes, max_genes)
  } else {
    subsample_idx <- NULL
  }
  
  # Initialize
  data_curr <- data
  t_opt <- NULL
  error_vec <- numeric(t_max)
  
  for (i in seq_len(t_max)) {
    data_prev <- data_curr
    data_curr <- as.matrix(diff_op %*% data_curr)
    
    # Calculate Procrustes disparity
    if (!is.null(subsample_idx)) {
      error <- procrustes_disparity(
        data_prev[, subsample_idx, drop = FALSE],
        data_curr[, subsample_idx, drop = FALSE]
      )
    } else {
      error <- procrustes_disparity(data_prev, data_curr)
    }
    
    error_vec[i] <- error
    
    # Check convergence
    if (is.null(t_opt) && error < threshold) {
      t_opt <- i
      vmessage(sprintf("  Optimal t: %d (disparity=%.6f)", t_opt, error), verbose)
      break
    }
  }
  
  # Use t_max if no convergence
  if (is.null(t_opt)) {
    t_opt <- t_max
    warning(sprintf("Optimal t not found within t_max=%d. Using t_max.", t_max))
  }
  
  # Preserve names
  rownames(data_curr) <- rownames(data)
  colnames(data_curr) <- colnames(data)
  
  # Store metadata
  attr(data_curr, "optimal_t") <- t_opt
  attr(data_curr, "error_vec") <- error_vec[seq_len(min(t_opt + 5, t_max))]
  
  return(data_curr)
}

#' Compute Optimal t
#'
#' @description Determines optimal diffusion steps using Procrustes analysis.
#'
#' @param data Matrix (cells x genes).
#' @param diff_op Diffusion operator.
#' @param t_max Integer. Maximum t to test. Default is 20.
#' @param threshold Numeric. Convergence threshold. Default is 0.001.
#' @param subsample_genes Integer or NULL. Genes to subsample. Default is 500.
#' @param plot Logical. Create diagnostic plot. Default is FALSE.
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A list with t_opt, error_vec, and threshold.
#'
#' @export
magic_optimal_t <- function(data,
                            diff_op,
                            t_max = 20,
                            threshold = 0.001,
                            subsample_genes = 500,
                            plot = FALSE,
                            verbose = FALSE) {
  
  check_positive_int(t_max, "t_max")
  check_positive(threshold, "threshold")
  
  n_genes <- ncol(data)
  
  # Subsample genes
  if (!is.null(subsample_genes) && n_genes > subsample_genes) {
    gene_idx <- sample(n_genes, subsample_genes)
    data_sub <- data[, gene_idx, drop = FALSE]
  } else {
    data_sub <- data
  }
  
  vmessage(sprintf("Computing optimal t (t_max=%d, threshold=%.4f)...",
                   t_max, threshold), verbose)
  
  # Initialize
  data_curr <- to_dense(data_sub)
  error_vec <- numeric(t_max)
  t_opt <- NULL
  
  for (i in seq_len(t_max)) {
    data_prev <- data_curr
    data_curr <- as.matrix(diff_op %*% data_curr)
    error_vec[i] <- procrustes_disparity(data_prev, data_curr)
    
    if (is.null(t_opt) && error_vec[i] < threshold) {
      t_opt <- i
    }
  }
  
  if (is.null(t_opt)) {
    t_opt <- t_max
    warning(sprintf("Optimal t not found. Using t_max=%d.", t_max))
  }
  
  vmessage(sprintf("  Optimal t: %d", t_opt), verbose)
  
  # Plot if requested
  if (plot) {
    plot_disparity(error_vec, t_opt, threshold)
  }
  
  return(list(
    t_opt = t_opt,
    error_vec = error_vec,
    threshold = threshold
  ))
}

#' Plot disparity curve
#' @param error_vec Vector of disparities
#' @param t_opt Optimal t
#' @param threshold Threshold value
#' @keywords internal
plot_disparity <- function(error_vec, t_opt, threshold) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(t = seq_along(error_vec), disparity = error_vec)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = disparity)) +
      ggplot2::geom_line(color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(color = "steelblue", size = 2) +
      ggplot2::geom_point(data = df[t_opt, ], color = "red", size = 4) +
      ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = "MAGIC: Optimal t Selection",
        subtitle = sprintf("Optimal t = %d", t_opt),
        x = "Diffusion steps (t)",
        y = "Procrustes Disparity"
      ) +
      ggplot2::theme_minimal()
    
    print(p)
    invisible(p)
  } else {
    # Base R plot
    plot(seq_along(error_vec), error_vec, type = "b", pch = 19,
         xlab = "t", ylab = "Procrustes Disparity",
         main = "MAGIC: Optimal t Selection")
    points(t_opt, error_vec[t_opt], col = "red", pch = 19, cex = 2)
    abline(h = threshold, lty = 2, col = "gray")
  }
}

#' Matrix Power
#'
#' @description Computes power of a matrix using binary exponentiation.
#'
#' @param M Square matrix
#' @param n Power (positive integer)
#'
#' @return M^n
#' @keywords internal
matrix_power <- function(M, n) {
  if (n == 0) {
    return(Matrix::Diagonal(nrow(M)))
  }
  if (n == 1) {
    return(M)
  }
  
  # Binary exponentiation: O(log n) multiplications
  result <- Matrix::Diagonal(nrow(M))
  base <- M
  
  while (n > 0) {
    if (n %% 2 == 1) {
      result <- result %*% base
    }
    base <- base %*% base
    n <- n %/% 2
  }
  
  return(result)
}
