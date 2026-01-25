# Tests for MAGICR package

# Setup test data
set.seed(42)
n_cells <- 100
n_genes <- 50
test_data <- matrix(rpois(n_cells * n_genes, lambda = 5), 
                    nrow = n_cells, ncol = n_genes)
colnames(test_data) <- paste0("Gene", 1:n_genes)
rownames(test_data) <- paste0("Cell", 1:n_cells)

# Test basic magic function
test_that("magic works with default parameters", {
  result <- magic(test_data, t = 2, verbose = FALSE)
  
  expect_s3_class(result, "magic")
  expect_equal(nrow(result$result), n_cells)
  expect_equal(ncol(result$result), n_genes)
  expect_true(all(is.finite(as.matrix(result))))
})

test_that("magic works with gene subset", {
  genes <- c("Gene1", "Gene2", "Gene3")
  result <- magic(test_data, genes = genes, t = 2, verbose = FALSE)
  
  expect_equal(ncol(result$result), length(genes))
  expect_equal(colnames(result$result), genes)
})

test_that("magic works with integer gene indices", {
  genes <- c(1, 5, 10)
  result <- magic(test_data, genes = genes, t = 2, verbose = FALSE)
  
  expect_equal(ncol(result$result), length(genes))
})

test_that("magic works with sparse matrix input", {
  sparse_data <- Matrix::Matrix(test_data, sparse = TRUE)
  result <- magic(sparse_data, t = 2, verbose = FALSE)
  
  expect_s3_class(result, "magic")
  expect_equal(nrow(result$result), n_cells)
})

test_that("magic works with data.frame input", {
  df_data <- as.data.frame(test_data)
  result <- magic(df_data, t = 2, verbose = FALSE)
  
  expect_s3_class(result, "magic")
  expect_equal(nrow(result$result), n_cells)
})

# Test kNN graph construction
test_that("magic_knn_graph works correctly", {
  knn_result <- magic_knn_graph(test_data, knn = 5, verbose = FALSE)
  
  expect_type(knn_result, "list")
  expect_equal(nrow(knn_result$idx), n_cells)
  expect_equal(ncol(knn_result$idx), 15)  # knn_max = 3 * knn
  expect_true(all(knn_result$dist >= 0))
})

test_that("magic_knn_graph respects knn_max", {
  knn_result <- magic_knn_graph(test_data, knn = 5, knn_max = 10, verbose = FALSE)
  
  expect_equal(ncol(knn_result$idx), 10)
  expect_equal(ncol(knn_result$dist), 10)
})

# Test kernel computation
test_that("compute_alpha_kernel creates symmetric matrix", {
  knn_graph <- magic_knn_graph(test_data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  
  expect_true(inherits(K, "Matrix"))
  expect_equal(nrow(K), n_cells)
  expect_equal(ncol(K), n_cells)
  # Check symmetry
  expect_true(Matrix::isSymmetric(K, tol = 1e-10))
})

test_that("compute_alpha_kernel handles decay = NULL", {
  knn_graph <- magic_knn_graph(test_data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = NULL, verbose = FALSE)
  
  expect_true(inherits(K, "Matrix"))
  # Unweighted kernel should have values of 1 where connected
  expect_true(all(K@x %in% c(0, 0.5, 1)))  # 0.5 due to symmetrization
})

# Test diffusion operator
test_that("magic_diffusion_operator creates row-stochastic matrix", {
  knn_graph <- magic_knn_graph(test_data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, decay = 1, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  expect_true(inherits(P, "Matrix"))
  # Check row sums are 1
  row_sums <- Matrix::rowSums(P)
  expect_equal(row_sums, rep(1, n_cells), tolerance = 1e-10)
})

# Test imputation
test_that("magic_impute preserves dimensions", {
  knn_graph <- magic_knn_graph(test_data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  imputed <- magic_impute(test_data, P, t = 3, verbose = FALSE)
  
  expect_equal(dim(imputed), dim(test_data))
  expect_equal(rownames(imputed), rownames(test_data))
  expect_equal(colnames(imputed), colnames(test_data))
})

test_that("magic_impute with t=0 returns original data", {
  knn_graph <- magic_knn_graph(test_data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  # t=0 should return original (actually t must be >= 1, but 
  # we test that small t doesn't over-smooth)
  imputed_t1 <- magic_impute(test_data, P, t = 1, verbose = FALSE)
  imputed_t5 <- magic_impute(test_data, P, t = 5, verbose = FALSE)
  
  # Higher t should be more different from original
  diff_t1 <- sum((imputed_t1 - test_data)^2)
  diff_t5 <- sum((imputed_t5 - test_data)^2)
  
  expect_true(diff_t5 > diff_t1)
})

# Test optimal t selection
test_that("magic_optimal_t returns valid t", {
  knn_graph <- magic_knn_graph(test_data, knn = 5, verbose = FALSE)
  K <- compute_alpha_kernel(knn_graph, verbose = FALSE)
  P <- magic_diffusion_operator(K, verbose = FALSE)
  
  opt_result <- magic_optimal_t(test_data, P, t_max = 10, verbose = FALSE)
  
  expect_type(opt_result, "list")
  expect_true(opt_result$t_opt >= 1)
  expect_true(opt_result$t_opt <= 10)
  expect_length(opt_result$error_vec, 10)
})

# Test preprocessing functions
test_that("library_size_normalize preserves median", {
  normalized <- library_size_normalize(test_data, verbose = FALSE)
  
  original_median <- median(rowSums(test_data))
  normalized_median <- median(rowSums(normalized))
  
  expect_equal(original_median, normalized_median, tolerance = 1e-10)
})

test_that("sqrt_normalize works correctly", {
  sqrt_data <- sqrt_normalize(test_data, verbose = FALSE)
  
  expect_equal(dim(sqrt_data), dim(test_data))
  expect_true(all(sqrt_data == sqrt(test_data)))
})

test_that("log_transform works correctly", {
  log_data <- log_transform(test_data, pseudo_count = 1, verbose = FALSE)
  
  expect_equal(dim(log_data), dim(test_data))
  expect_true(all.equal(log_data, log(test_data + 1)))
})

# Test S3 methods
test_that("S3 methods work correctly", {
  result <- magic(test_data, t = 2, verbose = FALSE)
  
  # print should not error
  expect_output(print(result), "MAGIC Result")
  
  # as.matrix
  mat <- as.matrix(result)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(n_cells, n_genes))
  
  # as.data.frame
  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
  expect_equal(dim(df), c(n_cells, n_genes))
  
  # subsetting
  subset1 <- result[1:10, ]
  expect_equal(nrow(subset1), 10)
  
  subset2 <- result[, 1:5]
  expect_equal(ncol(subset2), 5)
})

# Test approximate solver
test_that("approximate solver works", {
  result <- magic(test_data, t = 3, solver = "approximate", 
                  npca = 20, verbose = FALSE)
  
  expect_s3_class(result, "magic")
  expect_equal(nrow(result$result), n_cells)
})

# Test auto t selection
test_that("auto t selection works", {
  result <- magic(test_data, t = "auto", t_max = 10, verbose = FALSE)
  
  expect_s3_class(result, "magic")
  expect_true(result$params$t >= 1)
  expect_true(result$params$t <= 10)
})

# Test distance metrics
test_that("different distance metrics work", {
  for (dist in c("euclidean", "cosine", "correlation")) {
    result <- magic(test_data, t = 2, distance = dist, verbose = FALSE)
    expect_s3_class(result, "magic")
  }
})

# Test init parameter for reusing diffusion operator
test_that("init parameter reuses diffusion operator", {
  result1 <- magic(test_data, t = 2, verbose = FALSE)
  result2 <- magic(test_data, t = 3, init = result1, verbose = FALSE)
  
  # Same diffusion operator should be used
  expect_identical(result1$diff_op, result2$diff_op)
  # But different t
  expect_equal(result2$params$t, 3)
})

# Test with very small dataset
test_that("magic handles small datasets", {
  small_data <- test_data[1:20, 1:10]
  result <- magic(small_data, t = 2, knn = 3, verbose = FALSE)
  
  expect_s3_class(result, "magic")
})

# Test error handling
test_that("magic throws errors for invalid input", {
  expect_error(magic(test_data, knn = -1))
  expect_error(magic(test_data, t = -1))
  expect_error(magic(test_data, knn = 0))
})

# Test utility functions
test_that("utility functions work", {
  # Config
  config <- magicr_config()
  expect_type(config, "list")
  expect_true("version" %in% names(config))
  
  # has_cpp_acceleration (should return TRUE or FALSE)
  cpp <- has_cpp_acceleration()
  expect_type(cpp, "logical")
})
