#' @title Kernel and Graph Construction
#' @description Functions for constructing the kNN graph and diffusion kernel.
#' @name kernel
NULL

#' Build k-Nearest Neighbors Graph
#'
#' @description Constructs a k-nearest neighbors graph from the input data,
#'   computing distances and neighbor indices efficiently.
#'
#' @param data A matrix (cells x features) of expression values or PCA coordinates.
#' @param knn Integer. Number of nearest neighbors for adaptive bandwidth.
#'   Default is 5.
#' @param knn_max Integer or NULL. Maximum number of neighbors to consider.
#'   If NULL, defaults to 3 * knn.
#' @param distance Character. Distance metric. One of "euclidean" (default),
#'   "cosine", "manhattan", or "correlation".
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A list containing:
#'   \item{idx}{Matrix of neighbor indices (n_cells x knn_max)}
#'   \item{dist}{Matrix of distances to neighbors (n_cells x knn_max)}
#'   \item{knn}{Number of neighbors for bandwidth}
#'   \item{knn_max}{Maximum number of neighbors}
#'
#' @export
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(42)
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#'
#' # Build kNN graph
#' knn_result <- magic_knn_graph(data, knn = 5)
#' }
magic_knn_graph <- function(data,
                            knn = 5,
                            knn_max = NULL,
                            distance = c("euclidean", "cosine", "manhattan", "correlation"),
                            verbose = FALSE) {
  
  # Validate parameters
  check_positive_int(knn, "knn")
  distance <- match.arg(distance)
  
  if (is.null(knn_max)) {
    knn_max <- knn * 3
  }
  check_positive_int(knn_max, "knn_max")
  
  if (knn_max < knn) {
    knn_max <- knn
    vmessage("  knn_max adjusted to equal knn", verbose)
  }
  
  n_cells <- nrow(data)
  
  # Ensure knn_max doesn't exceed number of cells
  if (knn_max >= n_cells) {
    knn_max <- n_cells - 1
    vmessage(sprintf("  knn_max adjusted to %d (n_cells - 1)", knn_max), verbose)
  }
  
  # Convert to dense matrix for distance calculations
  data <- to_dense(data)
  
  vmessage(sprintf("Building kNN graph (k=%d, knn_max=%d, distance=%s)...",
                   knn, knn_max, distance), verbose)
  
  # Preprocess data based on distance metric
  if (distance == "cosine") {
    # Normalize rows to unit length for cosine distance
    row_norms <- sqrt(rowSums(data^2))
    row_norms[row_norms == 0] <- 1  # Avoid division by zero
    data <- data / row_norms
    # Now Euclidean distance on normalized data = sqrt(2 * (1 - cosine_similarity))
    use_distance <- "euclidean"
  } else if (distance == "correlation") {
    # Center and normalize each row for correlation distance
    data <- t(scale(t(data), center = TRUE, scale = TRUE))
    data[!is.finite(data)] <- 0
    use_distance <- "euclidean"
  } else {
    use_distance <- distance
  }
  
  # Use RANN for fast kNN search
  # RANN uses Euclidean distance by default, Manhattan also supported
  if (use_distance == "manhattan") {
    # RANN doesn't support Manhattan directly, use a workaround
    # or fall back to slower method
    vmessage("  Using Manhattan distance (may be slower)", verbose)
    knn_result <- knn_search_manhattan(data, knn_max + 1)
  } else {
    # Euclidean distance (includes transformed cosine and correlation)
    knn_result <- RANN::nn2(data, data, k = knn_max + 1, 
                            treetype = "kd", searchtype = "standard")
  }
  
  # Remove self (first neighbor is always the point itself)
  idx <- knn_result$nn.idx[, -1, drop = FALSE]
  dist <- knn_result$nn.dists[, -1, drop = FALSE]
  
  vmessage(sprintf("  kNN graph constructed: %d cells, %d neighbors each",
                   n_cells, ncol(idx)), verbose)
  
  return(list(
    idx = idx,
    dist = dist,
    knn = knn,
    knn_max = knn_max,
    distance = distance
  ))
}

#' Manhattan Distance kNN Search
#' @description Fallback kNN search using Manhattan distance.
#' @param data Matrix of data points
#' @param k Number of neighbors
#' @return List with nn.idx and nn.dists
#' @keywords internal
knn_search_manhattan <- function(data, k) {
  n <- nrow(data)
  k <- min(k, n)
  
  # Chunk processing for memory efficiency
  chunk_size <- min(1000, n)
  n_chunks <- ceiling(n / chunk_size)
  
  idx_list <- vector("list", n_chunks)
  dist_list <- vector("list", n_chunks)
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n)
    chunk_idx <- start_idx:end_idx
    
    # Compute Manhattan distances for this chunk
    chunk_data <- data[chunk_idx, , drop = FALSE]
    
    # Vectorized Manhattan distance computation
    # Using outer product approach in chunks to save memory
    chunk_dists <- matrix(0, nrow = length(chunk_idx), ncol = n)
    
    for (j in seq_len(ncol(data))) {
      chunk_dists <- chunk_dists + abs(outer(chunk_data[, j], data[, j], "-"))
    }
    
    # Find k nearest neighbors for each row
    idx_chunk <- t(apply(chunk_dists, 1, function(x) order(x)[1:k]))
    dist_chunk <- t(apply(seq_len(nrow(chunk_dists)), 1, function(row) {
      chunk_dists[row, idx_chunk[row, ]]
    }))
    
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
#' @description Computes the alpha-decaying kernel matrix from kNN distances.
#'   This is the core kernel used in MAGIC for diffusion.
#'
#' @param knn_graph A list from \code{magic_knn_graph} containing neighbor
#'   indices and distances.
#' @param decay Numeric. Alpha parameter for kernel decay rate. Higher values
#'   create sharper kernels. If NULL, uses a simple kNN kernel. Default is 1.
#' @param threshold Numeric. Threshold below which kernel values are set to 0.
#'   Default is 1e-4.
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A sparse symmetric kernel matrix (dgCMatrix).
#'
#' @details The alpha-decaying kernel is defined as:
#'   \deqn{K(i,j) = \exp\left(-\left(\frac{d(i,j)}{\epsilon_i}\right)^\alpha\right)}
#'   where \eqn{\epsilon_i} is the adaptive bandwidth (distance to the k-th neighbor)
#'   and \eqn{\alpha} is the decay parameter.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create example data and kNN graph
#' set.seed(42)
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' knn_graph <- magic_knn_graph(data, knn = 5)
#'
#' # Compute kernel
#' K <- compute_alpha_kernel(knn_graph, decay = 1)
#' }
compute_alpha_kernel <- function(knn_graph,
                                  decay = 1,
                                  threshold = 1e-4,
                                  verbose = FALSE) {
  
  idx <- knn_graph$idx
  dist <- knn_graph$dist
  knn <- knn_graph$knn
  
  n_cells <- nrow(idx)
  n_neighbors <- ncol(idx)
  
  vmessage(sprintf("Computing alpha-decaying kernel (decay=%.1f)...", 
                   decay %||% 0), verbose)
  
  if (is.null(decay) || decay == 0) {
    # Simple kNN kernel (unweighted)
    vmessage("  Using unweighted kNN kernel", verbose)
    K_values <- rep(1, n_cells * n_neighbors)
  } else {
    # Alpha-decaying kernel
    # Adaptive bandwidth: epsilon_i = distance to k-th neighbor
    epsilon <- dist[, knn]
    
    # Handle zero bandwidths
    epsilon[epsilon == 0] <- .Machine$double.eps
    
    # Compute kernel values: K(i,j) = exp(-(d(i,j)/epsilon_i)^decay)
    # Vectorized computation
    dist_normalized <- dist / epsilon  # Broadcast division
    K_values <- exp(-(dist_normalized^decay))
    
    # Apply threshold
    K_values[K_values < threshold] <- 0
    K_values <- as.vector(K_values)
  }
  
  # Build sparse matrix
  # Row indices
  i <- rep(seq_len(n_cells), times = n_neighbors)
  # Column indices
  j <- as.vector(idx)
  
  # Create sparse matrix
  K <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = K_values,
    dims = c(n_cells, n_cells)
  )
  
  # Symmetrize: K = (K + K') / 2
  # This makes the kernel undirected
  vmessage("  Symmetrizing kernel...", verbose)
  K <- (K + Matrix::t(K)) / 2
  
  # Report sparsity
  nnz <- Matrix::nnzero(K)
  density <- nnz / (n_cells^2)
  vmessage(sprintf("  Kernel matrix: %d x %d, %.2f%% non-zero",
                   n_cells, n_cells, density * 100), verbose)
  
  return(K)
}

#' Compute Diffusion Operator
#'
#' @description Constructs the row-stochastic diffusion operator (Markov matrix)
#'   from the kernel matrix.
#'
#' @param kernel A symmetric kernel matrix from \code{compute_alpha_kernel}.
#' @param verbose Logical or integer. Verbosity level.
#'
#' @return A row-stochastic sparse matrix (dgCMatrix) where each row sums to 1.
#'
#' @details The diffusion operator P is obtained by row-normalizing the kernel:
#'   \deqn{P = D^{-1} K}
#'   where D is the diagonal matrix of row sums of K.
#'
#'   This creates a Markov transition matrix representing the probability of
#'   transitioning between cells based on their similarity.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create kernel and diffusion operator
#' set.seed(42)
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' knn_graph <- magic_knn_graph(data, knn = 5)
#' K <- compute_alpha_kernel(knn_graph)
#' P <- magic_diffusion_operator(K)
#'
#' # Check row sums
#' all.equal(as.vector(Matrix::rowSums(P)), rep(1, nrow(P)))
#' }
magic_diffusion_operator <- function(kernel, verbose = FALSE) {
  vmessage("Computing diffusion operator (Markov normalization)...", verbose)
  
  # Row sums for normalization
  row_sums <- Matrix::rowSums(kernel)
  
  # Check for zero row sums
  if (any(row_sums == 0)) {
    warning("Some rows have zero sum (isolated points). Setting to uniform.")
    row_sums[row_sums == 0] <- 1
  }
  
  # Row normalization: P = D^(-1) * K
  # Use diagonal matrix multiplication for efficiency
  D_inv <- Matrix::Diagonal(x = 1 / row_sums)
  P <- D_inv %*% kernel
  
  # Ensure it's stored as dgCMatrix
  P <- as(P, "dgCMatrix")
  
  vmessage("  Diffusion operator computed", verbose)
  
  return(P)
}

#' Full Kernel Construction Pipeline
#'
#' @description Convenience function that builds the complete diffusion operator
#'   from raw data in one step.
#'
#' @param data Matrix of data (cells x features)
#' @param knn Number of neighbors for adaptive bandwidth
#' @param knn_max Maximum number of neighbors
#' @param decay Alpha decay parameter
#' @param distance Distance metric
#' @param threshold Kernel threshold
#' @param verbose Verbosity level
#'
#' @return A row-stochastic diffusion operator (sparse matrix)
#'
#' @keywords internal
build_diffusion_operator <- function(data,
                                      knn = 5,
                                      knn_max = NULL,
                                      decay = 1,
                                      distance = "euclidean",
                                      threshold = 1e-4,
                                      verbose = FALSE) {
  
  # Build kNN graph
  knn_graph <- magic_knn_graph(
    data = data,
    knn = knn,
    knn_max = knn_max,
    distance = distance,
    verbose = verbose
  )
  
  # Compute kernel
  kernel <- compute_alpha_kernel(
    knn_graph = knn_graph,
    decay = decay,
    threshold = threshold,
    verbose = verbose
  )
  
  # Build diffusion operator
  diff_op <- magic_diffusion_operator(
    kernel = kernel,
    verbose = verbose
  )
  
  return(diff_op)
}
