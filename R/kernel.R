#' @title Kernel and Graph Construction
#' @description Functions for kNN graph and diffusion kernel construction.
#' @name kernel
NULL

#' Build k-Nearest Neighbors Graph
#'
#' @description Constructs a kNN graph from input data.
#'
#' @param data Matrix (cells x features).
#' @param knn Integer. Number of neighbors for bandwidth. Default is 5.
#' @param knn_max Integer or NULL. Maximum neighbors. Default is 3 * knn.
#' @param distance Character. Distance metric: "euclidean", "cosine",
#'   "manhattan", or "correlation". Default is "euclidean".
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A list with idx (neighbor indices), dist (distances), and parameters.
#'
#' @export
magic_knn_graph <- function(data,
                            knn = 5,
                            knn_max = NULL,
                            distance = c("euclidean", "cosine", "manhattan", "correlation"),
                            verbose = FALSE) {
  
  check_positive_int(knn, "knn")
  distance <- match.arg(distance)
  
  if (is.null(knn_max)) {
    knn_max <- knn * 3
  }
  check_positive_int(knn_max, "knn_max")
  
  if (knn_max < knn) {
    knn_max <- knn
  }
  
  n_cells <- nrow(data)
  
  if (knn_max >= n_cells) {
    knn_max <- n_cells - 1
  }
  
  data <- to_dense(data)
  
  vmessage(sprintf("Building kNN graph (k=%d, knn_max=%d, distance=%s)...",
                   knn, knn_max, distance), verbose)
  
  # Transform data for different distance metrics
  if (distance == "cosine") {
    # L2 normalize rows
    row_norms <- sqrt(rowSums(data^2))
    row_norms[row_norms == 0] <- 1
    data <- data / row_norms
    use_distance <- "euclidean"
  } else if (distance == "correlation") {
    # Standardize rows
    data <- t(scale(t(data), center = TRUE, scale = TRUE))
    data[!is.finite(data)] <- 0
    use_distance <- "euclidean"
  } else {
    use_distance <- distance
  }
  
  # kNN search
  if (use_distance == "manhattan") {
    knn_result <- knn_manhattan(data, knn_max + 1)
  } else {
    knn_result <- RANN::nn2(data, data, k = knn_max + 1,
                            treetype = "kd", searchtype = "standard")
  }
  
  # Remove self (first neighbor)
  idx <- knn_result$nn.idx[, -1, drop = FALSE]
  dist <- knn_result$nn.dists[, -1, drop = FALSE]
  
  vmessage(sprintf("  Graph: %d cells, %d neighbors", n_cells, ncol(idx)), verbose)
  
  return(list(
    idx = idx,
    dist = dist,
    knn = knn,
    knn_max = knn_max,
    distance = distance
  ))
}

#' kNN with Manhattan distance
#' @param data Data matrix
#' @param k Number of neighbors
#' @return List with nn.idx and nn.dists
#' @keywords internal
knn_manhattan <- function(data, k) {
  n <- nrow(data)
  k <- min(k, n)
  
  chunk_size <- min(1000, n)
  n_chunks <- ceiling(n / chunk_size)
  
  idx_list <- vector("list", n_chunks)
  dist_list <- vector("list", n_chunks)
  
  for (i in seq_len(n_chunks)) {
    start <- (i - 1) * chunk_size + 1
    end <- min(i * chunk_size, n)
    chunk_idx <- start:end
    
    chunk_data <- data[chunk_idx, , drop = FALSE]
    chunk_dists <- matrix(0, nrow = length(chunk_idx), ncol = n)
    
    for (j in seq_len(ncol(data))) {
      chunk_dists <- chunk_dists + abs(outer(chunk_data[, j], data[, j], "-"))
    }
    
    idx_chunk <- t(apply(chunk_dists, 1, function(x) order(x)[seq_len(k)]))
    dist_chunk <- t(vapply(seq_len(nrow(chunk_dists)), function(row) {
      chunk_dists[row, idx_chunk[row, ]]
    }, numeric(k)))
    
    idx_list[[i]] <- idx_chunk
    dist_list[[i]] <- dist_chunk
  }
  
  list(
    nn.idx = do.call(rbind, idx_list),
    nn.dists = do.call(rbind, dist_list)
  )
}

#' Compute Alpha-Decaying Kernel
#'
#' @description Computes the alpha-decaying kernel from kNN distances.
#'
#' @param knn_graph Result from \code{magic_knn_graph}.
#' @param decay Numeric. Alpha parameter. If NULL, uses unweighted kernel. Default is 1.
#' @param threshold Numeric. Values below this are set to 0. Default is 1e-4.
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A sparse symmetric kernel matrix (dgCMatrix).
#'
#' @details The alpha-decaying kernel:
#'   \deqn{K(i,j) = \exp\left(-\left(\frac{d(i,j)}{\epsilon_i}\right)^\alpha\right)}
#'   
#'   where \eqn{\epsilon_i} is the distance to the k-th neighbor (adaptive bandwidth)
#'   and \eqn{\alpha} controls the decay rate.
#'   
#'   Reference: Moon, K.R., et al. (2019). Visualizing structure and transitions
#'   in high-dimensional biological data. Nature Biotechnology, 37(12), 1482-1492.
#'
#' @export
compute_alpha_kernel <- function(knn_graph,
                                 decay = 1,
                                 threshold = 1e-4,
                                 verbose = FALSE) {
  
  idx <- knn_graph$idx
  dist <- knn_graph$dist
  knn <- knn_graph$knn
  
  n_cells <- nrow(idx)
  n_neighbors <- ncol(idx)
  
  vmessage(sprintf("Computing kernel (decay=%.1f)...", decay %||% 0), verbose)
  
  if (is.null(decay) || decay == 0) {
    # Unweighted kNN kernel
    K_values <- rep(1, n_cells * n_neighbors)
  } else {
    # Alpha-decaying kernel
    # Bandwidth = distance to k-th neighbor
    epsilon <- dist[, knn]
    epsilon[epsilon == 0] <- .Machine$double.eps
    
    # K(i,j) = exp(-(d(i,j)/epsilon_i)^alpha)
    dist_scaled <- dist / epsilon
    K_values <- exp(-(dist_scaled^decay))
    K_values[K_values < threshold] <- 0
    K_values <- as.vector(K_values)
  }
  
  # Build sparse matrix
  i <- rep(seq_len(n_cells), times = n_neighbors)
  j <- as.vector(idx)
  
  K <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = K_values,
    dims = c(n_cells, n_cells)
  )
  
  # Symmetrize
  vmessage("  Symmetrizing kernel...", verbose)
  K <- (K + Matrix::t(K)) / 2
  
  nnz <- Matrix::nnzero(K)
  density <- nnz / (n_cells^2)
  vmessage(sprintf("  Kernel: %d x %d, %.2f%% non-zero", n_cells, n_cells, density * 100), verbose)
  
  return(K)
}

#' Compute Diffusion Operator
#'
#' @description Creates the row-stochastic diffusion operator (Markov matrix).
#'
#' @param kernel Symmetric kernel matrix from \code{compute_alpha_kernel}.
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A row-stochastic sparse matrix (each row sums to 1).
#'
#' @details The diffusion operator is:
#'   \deqn{P = D^{-1} K}
#'   where D is the diagonal matrix of row sums.
#'
#' @export
magic_diffusion_operator <- function(kernel, verbose = FALSE) {
  vmessage("Computing diffusion operator...", verbose)
  
  row_sums <- Matrix::rowSums(kernel)
  
  if (any(row_sums == 0)) {
    warning("Some cells have no neighbors. Setting to uniform distribution.")
    row_sums[row_sums == 0] <- 1
  }
  
  # P = D^{-1} * K
  D_inv <- Matrix::Diagonal(x = 1 / row_sums)
  P <- D_inv %*% kernel
  P <- as(P, "dgCMatrix")
  
  vmessage("  Diffusion operator computed", verbose)
  
  return(P)
}

#' Build Diffusion Operator (Pipeline)
#'
#' @description Builds diffusion operator from data in one step.
#'
#' @param data Matrix (cells x features)
#' @param knn Number of neighbors
#' @param knn_max Maximum neighbors
#' @param decay Alpha decay parameter
#' @param distance Distance metric
#' @param threshold Kernel threshold
#' @param verbose Verbosity level
#'
#' @return Row-stochastic diffusion operator
#'
#' @keywords internal
build_diffusion_operator <- function(data,
                                     knn = 5,
                                     knn_max = NULL,
                                     decay = 1,
                                     distance = "euclidean",
                                     threshold = 1e-4,
                                     verbose = FALSE) {
  
  knn_graph <- magic_knn_graph(
    data = data,
    knn = knn,
    knn_max = knn_max,
    distance = distance,
    verbose = verbose
  )
  
  kernel <- compute_alpha_kernel(
    knn_graph = knn_graph,
    decay = decay,
    threshold = threshold,
    verbose = verbose
  )
  
  diff_op <- magic_diffusion_operator(
    kernel = kernel,
    verbose = verbose
  )
  
  return(diff_op)
}
