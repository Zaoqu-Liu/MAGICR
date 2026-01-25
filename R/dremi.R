#' @title kNN-DREMI: Conditional Density Resampled Mutual Information
#' @description Functions for computing kNN-DREMI on single-cell data.
#' @name dremi
#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom graphics hist par
#' @importFrom stats var
NULL

#' Compute kNN-DREMI
#'
#' @description Calculates k-Nearest Neighbor conditional Density Resampled 
#'   Estimate of Mutual Information as defined in van Dijk et al, 2018.
#'
#' @details kNN-DREMI estimates the mutual information I(Y|X) using:
#'   1. Build kNN graph on X
#'   2. Resample Y based on local densities
#'   3. Compute mutual information on resampled data
#'   
#'   Note: kNN-DREMI is asymmetric - it estimates I(Y|X), not I(X|Y).
#'
#' @param x Numeric vector. Independent variable (gene X).
#' @param y Numeric vector. Dependent variable (gene Y).
#' @param k Integer. Number of nearest neighbors. Default is 10.
#' @param n_bins Integer. Number of bins for density resampling. Default is 20.
#' @param n_mesh Integer. Mesh points per bin. Default is 3.
#' @param n_jobs Integer. Parallel workers. Default is 1.
#' @param plot Logical. Create diagnostic plots. Default is FALSE.
#' @param return_drevi Logical. Also return DREVI values. Default is FALSE.
#'
#' @return If return_drevi is FALSE, returns the DREMI value (numeric).
#'   If return_drevi is TRUE, returns a list with dremi and drevi values.
#'
#' @references
#' van Dijk, D., et al. (2018). Recovering gene interactions from single-cell
#' data using data diffusion. Cell, 174(3), 716-729.
#' \doi{10.1016/j.cell.2018.05.061}
#'
#' @export
#' @examples
#' \dontrun{
#' # After running MAGIC
#' result <- magic(data, t = 3)
#' imputed <- as.matrix(result)
#' 
#' # Compute kNN-DREMI between two genes
#' dremi <- knnDREMI(imputed[, "GeneA"], imputed[, "GeneB"])
#' 
#' # With diagnostic plot
#' dremi <- knnDREMI(imputed[, "GeneA"], imputed[, "GeneB"], plot = TRUE)
#' }
knnDREMI <- function(x, y, 
                     k = 10, 
                     n_bins = 20, 
                     n_mesh = 3,
                     n_jobs = 1,
                     plot = FALSE,
                     return_drevi = FALSE) {
  
  # Validate inputs
  if (length(x) != length(y)) {
    stop("x and y must have the same length", call. = FALSE)
  }
  
  n <- length(x)
  if (n < k) {
    stop(sprintf("Not enough samples (%d) for k=%d neighbors", n, k), 
         call. = FALSE)
  }
  
  check_positive_int(k, "k")
  check_positive_int(n_bins, "n_bins")
  check_positive_int(n_mesh, "n_mesh")
  
  # Remove NAs
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]
  n <- length(x)
  
  # Normalize to [0, 1]
  x_norm <- (x - min(x)) / (max(x) - min(x) + .Machine$double.eps)
  y_norm <- (y - min(y)) / (max(y) - min(y) + .Machine$double.eps)
  
  # Build kNN on x
  x_mat <- matrix(x_norm, ncol = 1)
  knn_result <- RANN::nn2(x_mat, x_mat, k = k + 1)
  knn_idx <- knn_result$nn.idx[, -1, drop = FALSE]  # Remove self
  knn_dist <- knn_result$nn.dists[, -1, drop = FALSE]
  
  # Compute local density (bandwidth = k-th neighbor distance)
  bandwidth <- knn_dist[, k]
  bandwidth[bandwidth == 0] <- .Machine$double.eps
  
  # Create mesh grid
  x_mesh <- seq(0, 1, length.out = n_bins * n_mesh)
  y_mesh <- seq(0, 1, length.out = n_bins * n_mesh)
  
  # Compute conditional density P(Y|X) using kernel density estimation
  # For each x bin, estimate the density of y
  x_breaks <- seq(0, 1, length.out = n_bins + 1)
  x_mids <- (x_breaks[-1] + x_breaks[-(n_bins + 1)]) / 2
  
  # Initialize conditional density matrix
  cond_density <- matrix(0, nrow = n_bins, ncol = n_bins)
  
  for (i in seq_len(n_bins)) {
    # Find points in this x bin
    in_bin <- x_norm >= x_breaks[i] & x_norm < x_breaks[i + 1]
    if (i == n_bins) {
      in_bin <- x_norm >= x_breaks[i] & x_norm <= x_breaks[i + 1]
    }
    
    if (sum(in_bin) > 0) {
      y_in_bin <- y_norm[in_bin]
      
      # Histogram density for y in this x bin
      y_hist <- hist(y_in_bin, breaks = seq(0, 1, length.out = n_bins + 1), 
                     plot = FALSE)
      cond_density[i, ] <- y_hist$counts / sum(y_hist$counts + .Machine$double.eps)
    } else {
      # Uniform if no points
      cond_density[i, ] <- 1 / n_bins
    }
  }
  
  # Compute marginal density P(X)
  x_hist <- hist(x_norm, breaks = x_breaks, plot = FALSE)
  p_x <- x_hist$counts / sum(x_hist$counts + .Machine$double.eps)
  
  # Compute marginal density P(Y)
  y_breaks <- seq(0, 1, length.out = n_bins + 1)
  y_hist <- hist(y_norm, breaks = y_breaks, plot = FALSE)
  p_y <- y_hist$counts / sum(y_hist$counts + .Machine$double.eps)
  
  # Compute joint density P(X, Y)
  joint_density <- cond_density * p_x
  
  # Normalize joint density
  joint_density <- joint_density / (sum(joint_density) + .Machine$double.eps)
  
  # Compute mutual information
  # I(X;Y) = sum P(x,y) * log(P(x,y) / (P(x) * P(y)))
  mi <- 0
  for (i in seq_len(n_bins)) {
    for (j in seq_len(n_bins)) {
      if (joint_density[i, j] > 0 && p_x[i] > 0 && p_y[j] > 0) {
        mi <- mi + joint_density[i, j] * 
          log2(joint_density[i, j] / (p_x[i] * p_y[j] + .Machine$double.eps) + 
                 .Machine$double.eps)
      }
    }
  }
  
  # DREMI = rescaled MI to [0, 1]
  # Normalize by entropy of Y
  h_y <- -sum(p_y[p_y > 0] * log2(p_y[p_y > 0] + .Machine$double.eps))
  dremi <- mi / (h_y + .Machine$double.eps)
  dremi <- max(0, min(1, dremi))  # Clamp to [0, 1]
  
  # Compute DREVI if requested
  if (return_drevi) {
    # DREVI = variance explained
    # Compute E[Y|X] and Var[Y|X]
    y_mids <- (y_breaks[-1] + y_breaks[-(n_bins + 1)]) / 2
    
    e_y_given_x <- rowSums(cond_density * matrix(y_mids, nrow = n_bins, 
                                                  ncol = n_bins, byrow = TRUE))
    var_y_given_x <- rowSums(cond_density * 
                               (matrix(y_mids, nrow = n_bins, ncol = n_bins, 
                                       byrow = TRUE) - e_y_given_x)^2)
    
    # Expected conditional variance
    e_var_y_given_x <- sum(p_x * var_y_given_x)
    
    # Total variance of Y
    var_y <- var(y_norm)
    
    # DREVI = 1 - E[Var(Y|X)] / Var(Y)
    drevi <- 1 - e_var_y_given_x / (var_y + .Machine$double.eps)
    drevi <- max(0, min(1, drevi))
  }
  
  # Plot if requested
  if (plot) {
    plot_dremi(x, y, joint_density, dremi, 
               return_drevi && exists("drevi"), 
               if (return_drevi) drevi else NULL)
  }
  

  if (return_drevi) {
    return(list(dremi = dremi, drevi = drevi))
  } else {
    return(dremi)
  }
}

#' Plot DREMI diagnostic
#' @param x X values
#' @param y Y values  
#' @param joint_density Joint density matrix
#' @param dremi DREMI value
#' @param show_drevi Show DREVI
#' @param drevi DREVI value
#' @keywords internal
plot_dremi <- function(x, y, joint_density, dremi, show_drevi = FALSE, drevi = NULL) {
  # Avoid R CMD check note for ggplot aes variables
  density <- NULL
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Create data frame for scatter
    df <- data.frame(x = x, y = y)
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(alpha = 0.3, size = 0.5) +
      ggplot2::labs(
        title = sprintf("kNN-DREMI = %.3f%s", dremi,
                        if (show_drevi) sprintf(", DREVI = %.3f", drevi) else ""),
        x = "Gene X",
        y = "Gene Y"
      ) +
      ggplot2::theme_minimal()
    
    # Heatmap of joint density
    df_heat <- expand.grid(
      x = seq_len(nrow(joint_density)),
      y = seq_len(ncol(joint_density))
    )
    df_heat$density <- as.vector(joint_density)
    
    p2 <- ggplot2::ggplot(df_heat, ggplot2::aes(x = x, y = y, fill = density)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(
        title = "Joint Density P(X, Y)",
        x = "X bins",
        y = "Y bins"
      ) +
      ggplot2::theme_minimal()
    
    # Try to arrange side by side
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    } else {
      print(p1)
      print(p2)
    }
  } else {
    # Base R fallback
    par(mfrow = c(1, 2))
    
    plot(x, y, pch = 16, col = adjustcolor("steelblue", 0.3),
         main = sprintf("kNN-DREMI = %.3f", dremi),
         xlab = "Gene X", ylab = "Gene Y")
    
    image(joint_density, main = "Joint Density", 
          xlab = "X bins", ylab = "Y bins",
          col = hcl.colors(50, "viridis"))
    
    par(mfrow = c(1, 1))
  }
}

#' Compute kNN-DREMI on MAGIC result
#'
#' @description Convenience function to compute kNN-DREMI directly from
#'   a MAGIC result object.
#'
#' @param magic_result A magic object from \code{magic()}.
#' @param gene_x Name or index of gene X.
#' @param gene_y Name or index of gene Y.
#' @param ... Additional arguments passed to \code{knnDREMI}.
#'
#' @return DREMI value or list with DREMI and DREVI.
#'
#' @export
#' @examples
#' \dontrun{
#' result <- magic(data, t = 3)
#' dremi <- magic_knnDREMI(result, "GeneA", "GeneB")
#' }
magic_knnDREMI <- function(magic_result, gene_x, gene_y, ...) {
  if (!inherits(magic_result, "magic")) {
    stop("magic_result must be a magic object", call. = FALSE)
  }
  
  data <- as.matrix(magic_result)
  
  # Get gene indices/names
  if (is.character(gene_x)) {
    if (!gene_x %in% colnames(data)) {
      stop(sprintf("Gene '%s' not found", gene_x), call. = FALSE)
    }
    x <- data[, gene_x]
  } else {
    x <- data[, gene_x]
  }
  
  if (is.character(gene_y)) {
    if (!gene_y %in% colnames(data)) {
      stop(sprintf("Gene '%s' not found", gene_y), call. = FALSE)
    }
    y <- data[, gene_y]
  } else {
    y <- data[, gene_y]
  }
  
  knnDREMI(x, y, ...)
}
