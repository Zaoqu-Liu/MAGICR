#' @title kNN-DREMI: Conditional Density Resampled Mutual Information
#' @description Functions for computing kNN-DREMI on single-cell data.
#' @name dremi
#' @importFrom grDevices adjustcolor hcl.colors dev.copy dev.off png
#' @importFrom graphics hist par
#' @importFrom stats var sd
#' @importFrom rlang .data
NULL

#' Compute kNN-DREMI
#'
#' @description Calculates k-Nearest Neighbor conditional Density Resampled 
#'   Estimate of Mutual Information as defined in van Dijk et al, 2018.
#'
#' @details kNN-DREMI is an adaptation of DREMI (Krishnaswamy et al. 2014) for
#'   single cell RNA-sequencing data. DREMI captures the functional relationship
#'   between two genes across their entire dynamic range. The key change to
#'   kNN-DREMI is the replacement of the heat diffusion-based kernel-density
#'   estimator by a k-nearest neighbor-based density estimator.
#'   
#'   Note: kNN-DREMI, like Mutual Information and DREMI, is not symmetric.
#'   Here we are estimating I(Y|X).
#'   
#'   This implementation follows exactly the scprep Python implementation.
#'
#' @param x Numeric vector. Independent variable (gene X).
#' @param y Numeric vector. Dependent variable (gene Y).
#' @param k Integer. Number of nearest neighbors. Default is 10.
#' @param n_bins Integer. Number of bins for density resampling. Default is 20.
#' @param n_mesh Integer. Mesh density within bins. Default is 3.
#' @param n_jobs Integer. Number of parallel workers (not used in R). Default is 1.
#' @param plot Logical. Create diagnostic plots. Default is FALSE.
#' @param return_drevi Logical. Also return DREVI density matrix. Default is FALSE.
#' @param ... Additional arguments passed to plot function.
#'
#' @return If return_drevi is FALSE, returns the DREMI value (numeric).
#'   If return_drevi is TRUE, returns a list with dremi and drevi (density matrix).
#'
#' @references
#' van Dijk, D., et al. (2018). Recovering gene interactions from single-cell
#' data using data diffusion. Cell, 174(3), 716-729.
#' \doi{10.1016/j.cell.2018.05.061}
#' 
#' Krishnaswamy, S., et al. (2014). Conditional density-based analysis of T cell
#' signaling in single-cell data. Science, 346(6213), 1250689.
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
                     return_drevi = FALSE,
                     ...) {
  
  # Validate inputs
  x <- as.numeric(x)
  y <- as.numeric(y)
  
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
  
  # Check for constant arrays
  if (length(unique(x)) == 1 || length(unique(y)) == 1) {
    warning("Attempting to calculate kNN-DREMI on a constant array. Returning 0",
            call. = FALSE)
    if (return_drevi) {
      return(list(dremi = 0, drevi = matrix(0, n_bins, n_bins)))
    } else {
      return(0)
    }
  }
  
  # Step 0: Z-score normalize X and Y (exactly like Python)
  x_sd <- sd(x)
  y_sd <- sd(y)
  if (x_sd == 0) x_sd <- 1
  if (y_sd == 0) y_sd <- 1
  x <- (x - mean(x)) / x_sd
  y <- (y - mean(y)) / y_sd
  
  # Step 1: Create bin and mesh points
  # Note: n_bins + 1 edges for n_bins bins
  x_bins <- seq(min(x), max(x), length.out = n_bins + 1)
  y_bins <- seq(min(y), max(y), length.out = n_bins + 1)
  
  # Mesh points: ((n_mesh + 1) * n_bins) + 1 points
  n_mesh_points <- ((n_mesh + 1) * n_bins) + 1
  x_mesh <- seq(min(x), max(x), length.out = n_mesh_points)
  y_mesh <- seq(min(y), max(y), length.out = n_mesh_points)
  
  # Create mesh grid
  mesh_points <- expand.grid(x = x_mesh, y = y_mesh)
  mesh_points <- as.matrix(mesh_points)
  
  # Step 2: Find k nearest neighbors from mesh points to data
  data_xy <- cbind(x, y)
  knn_result <- RANN::nn2(data_xy, mesh_points, k = k)
  dists <- knn_result$nn.dists
  
  # Step 3: Compute density at each mesh point
  # density = k / (pi * r^2) where r is distance to k-th neighbor
  area <- pi * (dists[, k])^2
  area[area == 0] <- .Machine$double.eps
  density <- k / area
  
  # Step 4: Create mask for mesh points NOT on bin boundaries
  # (we exclude points exactly on bin edges)
  on_x_boundary <- mesh_points[, 1] %in% x_bins
  on_y_boundary <- mesh_points[, 2] %in% y_bins
  mesh_mask <- on_x_boundary | on_y_boundary
  
  # Step 5: Sum densities into bins using histogram2d
  # Only use mesh points that are NOT on boundaries
  valid_mesh <- !mesh_mask
  
  # Use hist2d equivalent
  bin_density <- matrix(0, n_bins, n_bins)
  
  if (sum(valid_mesh) > 0) {
    valid_x <- mesh_points[valid_mesh, 1]
    valid_y <- mesh_points[valid_mesh, 2]
    valid_density <- density[valid_mesh]
    
    # Assign each valid mesh point to its bin
    x_bin_idx <- findInterval(valid_x, x_bins, rightmost.closed = TRUE)
    y_bin_idx <- findInterval(valid_y, y_bins, rightmost.closed = TRUE)
    
    # Clamp indices
    x_bin_idx <- pmax(1, pmin(x_bin_idx, n_bins))
    y_bin_idx <- pmax(1, pmin(y_bin_idx, n_bins))
    
    # Sum densities into bins
    for (i in seq_along(valid_x)) {
      bin_density[x_bin_idx[i], y_bin_idx[i]] <- 
        bin_density[x_bin_idx[i], y_bin_idx[i]] + valid_density[i]
    }
  }
  
  # Transpose to match Python (histogram2d returns [x_bins, y_bins], we need [y_bins, x_bins])
  bin_density <- t(bin_density)
  
  # Normalize so total sums to 1
  total_density <- sum(bin_density)
  if (total_density > 0) {
    bin_density <- bin_density / total_density
  }
  
  # Step 6: Compute DREVI (conditional probability)
  # Normalize columns (each column sums to 1)
  col_sums <- colSums(bin_density)
  col_sums[col_sums == 0] <- 1  # Avoid division by zero
  drevi <- sweep(bin_density, 2, col_sums, "/")
  
  # Step 7: Compute conditional entropies for each column
  # H(Y|X=x) = -sum(p * log2(p))
  cond_entropies <- apply(drevi, 2, function(col) {
    col <- col[col > 0]
    if (length(col) == 0) return(0)
    -sum(col * log2(col))
  })
  
  # Step 8: Compute marginal entropy of Y
  p_y <- rowSums(bin_density)  # Marginal distribution of Y
  p_y <- p_y[p_y > 0]
  marginal_entropy <- -sum(p_y * log2(p_y))
  
  # Step 9: Compute mutual information (not normalized)
  # I(X;Y) = H(Y) - H(Y|X)
  # Conditional entropy = weighted average of column entropies
  p_x <- colSums(bin_density)  # Distribution of X (column sums)
  conditional_entropy <- sum(cond_entropies * p_x)
  mutual_info <- marginal_entropy - conditional_entropy
  
  # Step 10: Compute DREMI (normalized)
  # Use the normalized marginal entropy from DREVI
  # Note: scipy.stats.entropy normalizes the input, so we need to do the same
  p_y_norm <- rowSums(drevi)
  p_y_norm_sum <- sum(p_y_norm)
  if (p_y_norm_sum > 0) {
    p_y_norm <- p_y_norm / p_y_norm_sum  # Normalize to sum to 1
  }
  p_y_norm <- p_y_norm[p_y_norm > 0]
  marginal_entropy_norm <- -sum(p_y_norm * log2(p_y_norm))
  
  # Mean DREVI column for conditional entropy normalization
  # In Python: cond_sums_norm = np.mean(drevi) 
  cond_sums_norm <- mean(drevi)
  conditional_entropy_norm <- sum(cond_entropies * cond_sums_norm)
  
  dremi <- marginal_entropy_norm - conditional_entropy_norm
  
  # Clamp to reasonable range
  dremi <- max(0, dremi)
  
  # Plot if requested
  if (plot) {
    plot_knnDREMI(dremi, mutual_info, x, y, n_bins, n_mesh, 
                  density, bin_density, drevi, ...)
  }
  
  if (return_drevi) {
    return(list(dremi = dremi, drevi = drevi))
  } else {
    return(dremi)
  }
}

#' Plot kNN-DREMI diagnostic
#'
#' @description Creates diagnostic plots for kNN-DREMI like those seen in
#'   Fig 5C/D of van Dijk et al. 2018.
#'
#' @param dremi DREMI value
#' @param mutual_info Mutual information value
#' @param x Z-scored X values
#' @param y Z-scored Y values
#' @param n_bins Number of bins
#' @param n_mesh Mesh density
#' @param density Mesh point densities
#' @param bin_density Binned density matrix
#' @param drevi DREVI conditional probability matrix
#' @param figsize Figure size (ignored in R)
#' @param filename File to save plot (optional)
#' @param xlabel X-axis label
#' @param ylabel Y-axis label
#' @param ... Additional arguments
#'
#' @keywords internal
plot_knnDREMI <- function(dremi, mutual_info, x, y, n_bins, n_mesh,
                          density, bin_density, drevi,
                          figsize = c(12, 3.5),
                          filename = NULL,
                          xlabel = "Feature 1",
                          ylabel = "Feature 2",
                          ...) {
  
  if (requireNamespace("ggplot2", quietly = TRUE) && 
      requireNamespace("gridExtra", quietly = TRUE)) {
    
    # Plot 1: Raw data
    df1 <- data.frame(x = x, y = y)
    p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(size = 0.5, color = "black") +
      ggplot2::labs(title = "Input\ndata", x = xlabel, y = ylabel) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
    
    # Plot 2: kNN Density (mesh grid)
    n_mesh_points <- ((n_mesh + 1) * n_bins) + 1
    dens_mat <- matrix(log(density + 1e-10), nrow = n_mesh_points, ncol = n_mesh_points)
    df2 <- expand.grid(x = seq_len(n_mesh_points), y = seq_len(n_mesh_points))
    df2$density <- as.vector(dens_mat)
    
    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$density)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "inferno") +
      ggplot2::labs(title = "kNN\nDensity", x = xlabel, y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.position = "none"
      )
    
    # Plot 3: Joint probability
    df3 <- expand.grid(x = seq_len(n_bins), y = seq_len(n_bins))
    df3$density <- as.vector(bin_density)
    
    p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$density)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "inferno") +
      ggplot2::labs(
        title = sprintf("Joint Prob.\nMI=%.2f", mutual_info),
        x = xlabel, y = ""
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.position = "none"
      )
    
    # Plot 4: Conditional probability (DREVI)
    df4 <- expand.grid(x = seq_len(n_bins), y = seq_len(n_bins))
    df4$density <- as.vector(drevi)
    
    p4 <- ggplot2::ggplot(df4, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$density)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "inferno") +
      ggplot2::labs(
        title = sprintf("Conditional Prob.\nDREMI=%.2f", dremi),
        x = xlabel, y = ""
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.position = "none"
      )
    
    # Combine plots
    combined <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 4)
    
    if (!is.null(filename)) {
      ggplot2::ggsave(filename, combined, width = figsize[1], height = figsize[2])
    }
    
  } else {
    # Base R fallback
    old_par <- par(mfrow = c(1, 4), mar = c(3, 3, 3, 1))
    on.exit(par(old_par))
    
    # Plot 1: Raw data
    plot(x, y, pch = 16, cex = 0.3, col = "black",
         main = "Input\ndata", xlab = xlabel, ylab = ylabel,
         xaxt = "n", yaxt = "n")
    
    # Plot 2: kNN Density
    n_mesh_points <- ((n_mesh + 1) * n_bins) + 1
    dens_mat <- matrix(log(density + 1e-10), nrow = n_mesh_points, ncol = n_mesh_points)
    image(dens_mat, main = "kNN\nDensity", xlab = xlabel, ylab = "",
          col = hcl.colors(50, "Inferno"), xaxt = "n", yaxt = "n")
    
    # Plot 3: Joint probability
    image(bin_density, main = sprintf("Joint Prob.\nMI=%.2f", mutual_info),
          xlab = xlabel, ylab = "",
          col = hcl.colors(50, "Inferno"), xaxt = "n", yaxt = "n")
    
    # Plot 4: Conditional probability
    image(drevi, main = sprintf("Conditional Prob.\nDREMI=%.2f", dremi),
          xlab = xlabel, ylab = "",
          col = hcl.colors(50, "Inferno"), xaxt = "n", yaxt = "n")
    
    if (!is.null(filename)) {
      dev.copy(png, filename, width = figsize[1] * 100, height = figsize[2] * 100)
      dev.off()
    }
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
