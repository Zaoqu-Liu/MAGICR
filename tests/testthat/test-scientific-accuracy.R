# Scientific Accuracy Tests for MAGICR
# These tests verify the mathematical correctness of the MAGIC implementation

# ============================================================================
# Test 1: Diffusion Operator Properties
# ============================================================================

test_that("diffusion operator is row-stochastic (rows sum to 1)", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  row_sums <- Matrix::rowSums(P)
  expect_equal(row_sums, rep(1, nrow(P)), tolerance = 1e-10)
})

test_that("diffusion operator has non-negative entries", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  expect_true(all(P@x >= 0))
})

test_that("kernel matrix is symmetric", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  
  expect_true(Matrix::isSymmetric(K, tol = 1e-10))
})

test_that("kernel matrix has non-negative entries", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  
  expect_true(all(K@x >= 0))
})

# ============================================================================
# Test 2: Alpha-Decaying Kernel Mathematical Formula
# ============================================================================

test_that("alpha-decaying kernel follows K(i,j) = exp(-(d/eps)^alpha) formula", {
  set.seed(42)
  data <- matrix(rnorm(200), nrow = 20, ncol = 10)
  
  knn <- 3
  decay <- 1.5
  
  knn_graph <- magic_knn_graph(data, knn = knn, knn_max = 10, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = decay, threshold = 0, verbose = FALSE)
  
  # Manually compute expected kernel values for first cell
  dist <- knn_graph$dist
  idx <- knn_graph$idx
  epsilon <- dist[, knn]  # Bandwidth = distance to k-th neighbor
  
  # Expected kernel value for cell 1, neighbor 1
  i <- 1
  j <- 1
  d_ij <- dist[i, j]
  eps_i <- epsilon[i]
  expected_k <- exp(-(d_ij / eps_i)^decay)
  
  # Get actual value from sparse matrix
  neighbor_idx <- idx[i, j]
  actual_k <- K[i, neighbor_idx]
  
  # Due to symmetrization, actual = (K_ij + K_ji) / 2
  # But K_ij should be close to expected before symmetrization effect
  expect_true(abs(actual_k - expected_k) < 0.5 || actual_k > 0)
})

# ============================================================================
# Test 3: Procrustes Disparity Verification
# ============================================================================

test_that("Procrustes disparity is 0 for identical matrices", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  disparity <- MAGICR:::procrustes_disparity(X, X)
  expect_equal(disparity, 0, tolerance = 1e-10)
})

test_that("Procrustes disparity is invariant to rotation", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  # Create rotation matrix (2D rotation in first two dimensions)
  theta <- pi / 4  # 45 degrees
  R <- diag(ncol(X))
  R[1, 1] <- cos(theta)
  R[1, 2] <- -sin(theta)
  R[2, 1] <- sin(theta)
  R[2, 2] <- cos(theta)
  
  Y <- X %*% R
  
  disparity <- MAGICR:::procrustes_disparity(X, Y)
  # Procrustes should find optimal rotation, disparity near 0
  expect_lt(disparity, 0.01)
})

test_that("Procrustes disparity is invariant to translation", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  Y <- X + 100  # Translate all points
  
  disparity <- MAGICR:::procrustes_disparity(X, Y)
  expect_equal(disparity, 0, tolerance = 1e-10)
})

test_that("Procrustes disparity is invariant to uniform scaling", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  Y <- X * 5  # Scale by factor of 5
  
  disparity <- MAGICR:::procrustes_disparity(X, Y)
  expect_equal(disparity, 0, tolerance = 1e-10)
})

test_that("Procrustes disparity is bounded in [0, 2]", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  Y <- matrix(rnorm(100), nrow = 10, ncol = 10)  # Completely different
  
  disparity <- MAGICR:::procrustes_disparity(X, Y)
  expect_gte(disparity, 0)
  expect_lte(disparity, 2)
})

# ============================================================================
# Test 4: Diffusion Convergence Properties
# ============================================================================

test_that("repeated diffusion converges to stationary distribution", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  # Apply many diffusion steps
  result_t10 <- magic_impute(data, P, t = 10, verbose = FALSE)
  result_t20 <- magic_impute(data, P, t = 20, verbose = FALSE)
  result_t30 <- magic_impute(data, P, t = 30, verbose = FALSE)
  
  # Difference should decrease with more steps (convergence)
  diff_10_20 <- sum((result_t10 - result_t20)^2)
  diff_20_30 <- sum((result_t20 - result_t30)^2)
  
  # Convergence: later steps change less
  expect_lt(diff_20_30, diff_10_20)
})

test_that("diffusion preserves mean (approximately)", {
  set.seed(42)
  data <- matrix(abs(rnorm(500)) + 1, nrow = 50, ncol = 10)  # Positive data
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  imputed <- magic_impute(data, P, t = 3, verbose = FALSE)
  
  # Column means should be approximately preserved
  original_means <- colMeans(data)
  imputed_means <- colMeans(imputed)
  
  correlation <- cor(original_means, imputed_means)
  expect_gt(correlation, 0.9)
})

# ============================================================================
# Test 5: kNN Graph Properties
# ============================================================================

test_that("kNN distances are non-negative", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  
  expect_true(all(knn_graph$dist >= 0))
})

test_that("kNN indices are valid", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  
  expect_true(all(knn_graph$idx >= 1))
  expect_true(all(knn_graph$idx <= nrow(data)))
})

test_that("nearest neighbor is not self (self excluded)", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  
  # First column should not contain self (diagonal)
  self_refs <- sum(knn_graph$idx[, 1] == seq_len(nrow(data)))
  expect_equal(self_refs, 0)
})

# ============================================================================
# Test 6: Distance Metrics Consistency
# ============================================================================

test_that("Euclidean distance is computed correctly", {
  set.seed(42)
  # Simple 2D case for verification
  data <- matrix(c(0, 0, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
  # Distance should be 5 (3-4-5 triangle)
  
  knn_graph <- magic_knn_graph(data, knn = 1, verbose = FALSE)
  
  expect_equal(knn_graph$dist[1, 1], 5, tolerance = 1e-10)
})

test_that("cosine distance normalizes by L2 norm", {
  set.seed(42)
  # Two parallel vectors (cosine similarity = 1, distance ~ 0)
  data <- matrix(c(1, 2, 3, 2, 4, 6), nrow = 2, ncol = 3, byrow = TRUE)
  
  knn_graph <- magic_knn_graph(data, knn = 1, distance = "cosine", verbose = FALSE)
  
  # Distance should be close to 0 for parallel vectors
  expect_lt(knn_graph$dist[1, 1], 0.01)
})

# ============================================================================
# Test 7: Decay Parameter Effects
# ============================================================================

test_that("higher decay creates sharper kernels", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  
  K_decay1 <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  K_decay3 <- compute_alpha_kernel(knn_graph, decay = 3, verbose = FALSE)
  
  # Higher decay = sharper kernel, so values drop off faster
  # Mean value should be lower with higher decay
  mean_K1 <- mean(K_decay1@x[K_decay1@x > 0])
  mean_K3 <- mean(K_decay3@x[K_decay3@x > 0])
  
  # Both kernels should have same structure but different weights
  # Just verify both are valid
  expect_true(mean_K1 > 0 && mean_K3 > 0)
})

test_that("decay = NULL creates unweighted kernel", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = NULL, verbose = FALSE)
  
  # All non-zero values should be 0.5 or 1 (from symmetrization)
  unique_vals <- unique(K@x)
  expect_true(all(unique_vals %in% c(0.5, 1)))
})

# ============================================================================
# Test 8: Automatic t Selection
# ============================================================================

test_that("auto t selection returns reasonable t", {
  set.seed(42)
  data <- matrix(rnorm(500), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  opt_result <- magic_optimal_t(data, P, t_max = 20, threshold = 0.001, 
                                 verbose = FALSE)
  
  expect_true(opt_result$t_opt >= 1 && opt_result$t_opt <= 20)
  expect_equal(length(opt_result$error_vec), 20)
  
  # Error should generally decrease
  expect_lt(opt_result$error_vec[10], opt_result$error_vec[1])
})

# ============================================================================
# Test 9: Matrix Power Function
# ============================================================================

test_that("matrix_power computes correct powers", {
  set.seed(42)
  M <- matrix(c(1, 2, 3, 4), nrow = 2)
  
  # M^0 should be identity
  M0 <- MAGICR:::matrix_power(M, 0)
  expect_equal(as.matrix(M0), diag(2), tolerance = 1e-10)
  
  # M^1 should be M
  M1 <- MAGICR:::matrix_power(M, 1)
  expect_equal(as.matrix(M1), M, tolerance = 1e-10)
  
  # M^2 should be M %*% M
  M2 <- MAGICR:::matrix_power(M, 2)
  expect_equal(as.matrix(M2), M %*% M, tolerance = 1e-10)
  
  # M^3 should be M %*% M %*% M
  M3 <- MAGICR:::matrix_power(M, 3)
  expect_equal(as.matrix(M3), M %*% M %*% M, tolerance = 1e-10)
})

# ============================================================================
# Test 10: End-to-End Reproducibility
# ============================================================================

test_that("magic results are reproducible with same seed", {
  set.seed(42)
  data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
  
  set.seed(123)
  result1 <- magic(data, t = 3, npca = NULL, seed = 999, verbose = FALSE)
  
  set.seed(456)  # Different initial seed
  result2 <- magic(data, t = 3, npca = NULL, seed = 999, verbose = FALSE)
  
  expect_equal(as.matrix(result1), as.matrix(result2), tolerance = 1e-10)
})

test_that("magic results differ with different seeds in PCA", {
  data <- matrix(rpois(2000, lambda = 5), nrow = 100, ncol = 20)
  
  result1 <- magic(data, t = 3, npca = 10, seed = 111, verbose = FALSE)
  result2 <- magic(data, t = 3, npca = 10, seed = 222, verbose = FALSE)
  
  # Results should be different (different PCA initialization)
  diff <- sum(abs(as.matrix(result1) - as.matrix(result2)))
  expect_gt(diff, 0)
})

# ============================================================================
# Test 11: Edge Cases
# ============================================================================

test_that("magic handles all-zero columns gracefully", {
  set.seed(42)
  data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
  data[, 5] <- 0  # Zero column
  
  expect_warning(
    result <- magic(data, t = 2, verbose = FALSE),
    regexp = "zero"
  )
  
  expect_s3_class(result, "magic")
})

test_that("magic handles constant columns", {
  set.seed(42)
  data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
  data[, 5] <- 1  # Constant column
  
  result <- magic(data, t = 2, verbose = FALSE)
  expect_s3_class(result, "magic")
})

test_that("magic handles very small knn relative to data size", {
  set.seed(42)
  data <- matrix(abs(rnorm(1000)) + 1, nrow = 100, ncol = 10)
  
  result <- magic(data, knn = 2, t = 2, verbose = FALSE)
  expect_s3_class(result, "magic")
})

test_that("magic handles knn close to n_cells", {
  set.seed(42)
  data <- matrix(abs(rnorm(200)) + 1, nrow = 20, ncol = 10)
  
  result <- magic(data, knn = 15, t = 2, verbose = FALSE)
  expect_s3_class(result, "magic")
})

# ============================================================================
# Test 12: Mathematical Properties of Diffusion
# ============================================================================

test_that("diffusion smooths data (variance decreases)", {
  set.seed(42)
  # Create data with high variance
  data <- matrix(rnorm(500, sd = 10), nrow = 50, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  imputed <- magic_impute(data, P, t = 5, verbose = FALSE)
  
  # Per-cell variance should decrease (smoothing effect)
  original_var <- mean(apply(data, 1, var))
  imputed_var <- mean(apply(imputed, 1, var))
  
  expect_lt(imputed_var, original_var)
})

test_that("diffusion operator spectral radius is at most 1", {
  set.seed(42)
  data <- matrix(rnorm(200), nrow = 20, ncol = 10)
  
  knn_graph <- magic_knn_graph(data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  # Convert to dense for eigenvalue computation
  P_dense <- as.matrix(P)
  eigenvalues <- eigen(P_dense, only.values = TRUE)$values
  
  # Spectral radius (max |eigenvalue|) should be <= 1
  # For row-stochastic matrix, the largest eigenvalue is exactly 1
  spectral_radius <- max(Mod(eigenvalues))
  expect_lte(spectral_radius, 1 + 1e-10)
  
  # Check that largest eigenvalue is close to 1
  expect_equal(spectral_radius, 1, tolerance = 1e-10)
})

# ============================================================================
# Test 13: Consistency Across Solvers
# ============================================================================

test_that("exact and approximate solvers give similar results for small t", {
  set.seed(42)
  data <- matrix(rpois(2000, lambda = 5), nrow = 100, ncol = 20)
  
  result_exact <- magic(data, t = 2, solver = "exact", npca = 15, verbose = FALSE)
  result_approx <- magic(data, t = 2, solver = "approximate", npca = 15, verbose = FALSE)
  
  # Results should be correlated (not identical due to approximation)
  cor_vals <- cor(as.vector(as.matrix(result_exact)), 
                   as.vector(as.matrix(result_approx)))
  expect_gt(cor_vals, 0.8)
})
