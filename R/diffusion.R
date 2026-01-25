#' @title Diffusion Operations
#' @description Functions for performing diffusion-based imputation.
#' @name diffusion
NULL

#' Perform MAGIC Imputation
#'
#' @description Applies the diffusion operator to impute/denoise gene expression.
#'
#' @param data Matrix (cells x genes) to impute.
#' @param diff_op Diffusion operator (Markov matrix) from \code{magic_diffusion_operator}.
#' @param t Integer. Number of diffusion steps. Higher values = more smoothing.
#'   If "auto", optimal t is selected automatically.
#' @param t_max Integer. Maximum t to consider when t = "auto". Default is 20.
#' @param threshold Numeric. Convergence threshold for automatic t selection.
#'   Default is 0.001.
#' @param solver Character. Either "exact" or "approximate". Default is "exact".
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A matrix of imputed values with the same dimensions as input.
#'
#' @details MAGIC imputation is performed by repeatedly applying the diffusion
#'   operator to the data:
#'   \deqn{X_{imputed} = P^t \times X}
#'   where P is the diffusion operator and t is the number of steps.
#'
#'   When t = "auto", the algorithm iteratively applies P and stops when the
#'   Procrustes disparity between successive iterations falls below the threshold.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(42)
#' data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#'
#' # Build diffusion operator
#' knn_graph <- magic_knn_graph(data, knn = 5)
#' K <- compute_alpha_kernel(knn_graph)
#' P <- magic_diffusion_operator(K)
#'
#' # Impute with fixed t
#' data_imputed <- magic_impute(data, P, t = 3)
#'
#' # Impute with automatic t selection
#' data_imputed_auto <- magic_impute(data, P, t = "auto")
#' }
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
    stop("diff_op must be a matrix or Matrix object", call. = FALSE)
  }
  
  if (nrow(diff_op) != ncol(diff_op)) {
    stop("diff_op must be square", call. = FALSE)
  }
  
  if (nrow(diff_op) != nrow(data)) {
    stop("Number of rows in data must match dimensions of diff_op", call. = FALSE)
  }
  
  # Handle t parameter
  auto_t <- identical(t, "auto")
  if (!auto_t) {
    check_positive_int(t, "t")
  }
  check_positive_int(t_max, "t_max")
  check_positive(threshold, "threshold")
  
  # Convert data to dense matrix for imputation
  data <- to_dense(data)
  
  if (auto_t) {
    result <- magic_impute_auto_t(
      data = data,
      diff_op = diff_op,
      t_max = t_max,
      threshold = threshold,
      verbose = verbose
    )
  } else {
    result <- magic_impute_fixed_t(
      data = data,
      diff_op = diff_op,
      t = t,
      verbose = verbose
    )
  }
  
  return(result)
}

#' Impute with Fixed t
#'
#' @description Internal function for imputation with fixed t.
#' @param data Data matrix
#' @param diff_op Diffusion operator
#' @param t Number of steps
#' @param verbose Verbosity level
#' @return Imputed matrix
#' @keywords internal
magic_impute_fixed_t <- function(data, diff_op, t, verbose = FALSE) {
  n_cells <- nrow(data)
  n_genes <- ncol(data)
  
  vmessage(sprintf("Imputing with t=%d...", t), verbose)
  
  # Decide strategy based on matrix sizes
  # If diff_op is smaller than data (n_cells < n_genes), power the operator first
  # Otherwise, iterate multiplication
  
  if (n_cells < n_genes && t > 1) {
    # Classic MAGIC: (P^t) * X is more efficient
    vmessage(sprintf("  Using matrix power strategy (n_cells=%d < n_genes=%d)",
                     n_cells, n_genes), verbose, level = 2)
    
    # Power the diffusion operator
    P_t <- matrix_power(diff_op, t, verbose = verbose)
    
    # Single matrix multiplication
    result <- as.matrix(P_t %*% data)
    
  } else {
    # Fast MAGIC: iterative P * (P * ... * X)
    vmessage("  Using iterative multiplication strategy", verbose, level = 2)
    
    result <- data
    for (i in seq_len(t)) {
      result <- as.matrix(diff_op %*% result)
      vmessage(sprintf("    Step %d/%d", i, t), verbose, level = 2)
    }
  }
  
  # Preserve names
  rownames(result) <- rownames(data)
  colnames(result) <- colnames(data)
  
  vmessage("  Imputation complete", verbose)
  
  return(result)
}

#' Impute with Automatic t Selection
#'
#' @description Internal function for imputation with automatic t selection.
#' @param data Data matrix
#' @param diff_op Diffusion operator
#' @param t_max Maximum t to try
#' @param threshold Convergence threshold
#' @param verbose Verbosity level
#' @return Imputed matrix
#' @keywords internal
magic_impute_auto_t <- function(data, diff_op, t_max, threshold, verbose = FALSE) {
  n_cells <- nrow(data)
  n_genes <- ncol(data)
  
  vmessage("Selecting optimal t automatically...", verbose)
  
  # Subsample genes for faster error calculation if necessary
  max_genes_for_error <- 500
  if (n_genes > max_genes_for_error) {
    subsample_idx <- sample(n_genes, max_genes_for_error)
    vmessage(sprintf("  Subsampling %d genes for optimal t calculation",
                     max_genes_for_error), verbose, level = 2)
  } else {
    subsample_idx <- NULL
  }
  
  # Initialize
  data_curr <- data
  data_prev <- NULL
  t_opt <- NULL
  error_vec <- numeric(t_max)
  
  for (i in seq_len(t_max)) {
    # Apply one step of diffusion
    data_prev <- data_curr
    data_curr <- as.matrix(diff_op %*% data_curr)
    
    # Calculate error on subsampled genes
    if (!is.null(subsample_idx)) {
      error <- procrustes_disparity(
        data_prev[, subsample_idx, drop = FALSE],
        data_curr[, subsample_idx, drop = FALSE]
      )
    } else {
      error <- procrustes_disparity(data_prev, data_curr)
    }
    
    error_vec[i] <- error
    vmessage(sprintf("    t=%d: disparity=%.6f", i, error), verbose, level = 2)
    
    # Check convergence
    if (is.null(t_opt) && error < threshold) {
      t_opt <- i
      vmessage(sprintf("  Optimal t selected: %d (disparity=%.6f)", t_opt, error),
               verbose)
      # Continue for a bit to ensure stability, but keep this result
      break
    }
  }
  
  # If no convergence, use t_max
  if (is.null(t_opt)) {
    t_opt <- t_max
    warning(sprintf("Optimal t not found within t_max=%d. Using t_max.", t_max))
    vmessage(sprintf("  Using t=%d (max)", t_opt), verbose)
  }
  
  # Preserve names
  rownames(data_curr) <- rownames(data)
  colnames(data_curr) <- colnames(data)
  
  # Store optimal t as attribute
  attr(data_curr, "optimal_t") <- t_opt
  attr(data_curr, "error_vec") <- error_vec[1:min(t_opt + 5, t_max)]
  
  return(data_curr)
}

#' Compute Optimal t for MAGIC
#'
#' @description Determines the optimal number of diffusion steps using
#'   Procrustes disparity analysis.
#'
#' @param data Matrix (cells x genes) of expression values.
#' @param diff_op Diffusion operator.
#' @param t_max Integer. Maximum t to test. Default is 20.
#' @param threshold Numeric. Convergence threshold. Default is 0.001.
#' @param subsample_genes Integer or NULL. Number of genes to subsample for
#'   faster computation. If NULL, uses all genes. Default is 500.
#' @param plot Logical. Whether to create a diagnostic plot. Default is FALSE.
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A list containing:
#'   \item{t_opt}{Optimal t value}
#'   \item{error_vec}{Vector of Procrustes disparities for each t}
#'   \item{threshold}{The threshold used}
#'
#' @export
#' @examples
#' \dontrun{
#' # Analyze optimal t
#' opt_t_result <- magic_optimal_t(data, P, t_max = 20, plot = TRUE)
#' print(opt_t_result$t_opt)
#' }
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
  
  # Subsample genes if necessary
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
    
    vmessage(sprintf("  t=%d: disparity=%.6f", i, error_vec[i]), verbose, level = 2)
    
    if (is.null(t_opt) && error_vec[i] < threshold) {
      t_opt <- i
    }
  }
  
  if (is.null(t_opt)) {
    t_opt <- t_max
    warning(sprintf("Optimal t not found below threshold. Using t_max=%d.", t_max))
  }
  
  vmessage(sprintf("  Optimal t: %d", t_opt), verbose)
  
  # Create plot if requested
  if (plot) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_optimal_t(error_vec, t_opt, threshold)
    } else {
      # Base R plot
      plot(seq_len(t_max), error_vec, type = "b", pch = 19,
           xlab = "t", ylab = "Procrustes Disparity",
           main = "Optimal t Selection")
      points(t_opt, error_vec[t_opt], col = "red", pch = 19, cex = 2)
      abline(h = threshold, lty = 2, col = "gray")
      legend("topright", legend = c(sprintf("Optimal t = %d", t_opt), 
                                     sprintf("Threshold = %.4f", threshold)),
             col = c("red", "gray"), pch = c(19, NA), lty = c(NA, 2))
    }
  }
  
  return(list(
    t_opt = t_opt,
    error_vec = error_vec,
    threshold = threshold
  ))
}

#' Plot Optimal t Analysis
#' @param error_vec Vector of errors
#' @param t_opt Optimal t
#' @param threshold Threshold value
#' @keywords internal
plot_optimal_t <- function(error_vec, t_opt, threshold) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }
  
  df <- data.frame(
    t = seq_along(error_vec),
    disparity = error_vec
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = disparity)) +
    ggplot2::geom_line(color = "steelblue", size = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::geom_point(data = df[t_opt, ], color = "red", size = 4) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "gray50") +
    ggplot2::annotate("text", x = length(error_vec) * 0.8, y = threshold * 1.5,
                      label = sprintf("threshold = %.4f", threshold)) +
    ggplot2::labs(
      title = "MAGIC Optimal t Selection",
      subtitle = sprintf("Optimal t = %d", t_opt),
      x = "Diffusion steps (t)",
      y = "Procrustes Disparity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  
  print(p)
  invisible(p)
}

#' Matrix Power
#'
#' @description Efficiently computes the power of a sparse matrix.
#'
#' @param M Square matrix
#' @param n Power (positive integer)
#' @param verbose Verbosity level
#'
#' @return M^n
#' @keywords internal
matrix_power <- function(M, n, verbose = FALSE) {
  if (n == 0) {
    return(Matrix::Diagonal(nrow(M)))
  }
  if (n == 1) {
    return(M)
  }
  
  # Binary exponentiation for efficiency
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
