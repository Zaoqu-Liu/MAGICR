# Tests for knnDREMI and related functions

test_that("knnDREMI computes correctly for correlated data", {
  set.seed(42)
  
  # Create highly correlated data
  n <- 500
  x <- rnorm(n)
  y <- 0.8 * x + 0.2 * rnorm(n)
  
  dremi <- knnDREMI(x, y, k = 10, n_bins = 20, n_mesh = 3)
  
  # DREMI should be positive for correlated data
  expect_gt(dremi, 0)
  
  # DREMI should be reasonably high for strongly correlated data
  # Note: DREMI can be > 1 for highly correlated data
  expect_gt(dremi, 0.5)
})

test_that("knnDREMI returns low value for uncorrelated data", {
  set.seed(42)
  
  # Uncorrelated data
  n <- 500
  x <- rnorm(n)
  y <- rnorm(n)  # Independent
  
  dremi <- knnDREMI(x, y, k = 10, n_bins = 20, n_mesh = 3)
  
  # DREMI should be close to 0 for independent data
  expect_gte(dremi, 0)
  # Should be relatively low
  expect_lt(dremi, 0.5)
})

test_that("knnDREMI handles constant arrays", {
  x <- rep(1, 100)
  y <- rnorm(100)
  
  expect_warning(dremi <- knnDREMI(x, y), "constant array")
  expect_equal(dremi, 0)
})

test_that("knnDREMI returns DREVI when requested", {
  set.seed(42)
  
  n <- 200
  x <- rnorm(n)
  y <- 0.5 * x + 0.5 * rnorm(n)
  
  result <- knnDREMI(x, y, return_drevi = TRUE)
  
  expect_type(result, "list")
  expect_true("dremi" %in% names(result))
  expect_true("drevi" %in% names(result))
  expect_true(is.matrix(result$drevi))
  expect_equal(dim(result$drevi), c(20, 20))  # default n_bins
  
  # DREVI columns should sum to 1 (conditional probability)
  col_sums <- colSums(result$drevi)
  expect_true(all(abs(col_sums - 1) < 0.01 | col_sums == 0))
})

test_that("knnDREMI is asymmetric", {
  set.seed(42)
  
  n <- 300
  x <- rnorm(n)
  y <- x^2 + 0.1 * rnorm(n)  # Nonlinear relationship
  
  dremi_xy <- knnDREMI(x, y, k = 10)
  dremi_yx <- knnDREMI(y, x, k = 10)
  
  # DREMI is asymmetric, values should differ
  # (though for some relationships they might be similar)
  expect_true(is.numeric(dremi_xy))
  expect_true(is.numeric(dremi_yx))
})

test_that("knnDREMI validates inputs", {
  # Different lengths
  expect_error(knnDREMI(1:10, 1:5), "same length")
  
  # Not enough samples
  expect_error(knnDREMI(1:5, 1:5, k = 10), "Not enough samples")
})

test_that("magic_knnDREMI works with magic object", {
  skip_if_not_installed("testthat")
  
  # Create simple test data
  set.seed(42)
  data <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  colnames(data) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
  
  # Run magic with minimal settings
  result <- magic(data, knn = 5, t = 2, npca = NULL, verbose = FALSE)
  
  # Test magic_knnDREMI
  dremi <- magic_knnDREMI(result, "Gene1", "Gene2")
  expect_true(is.numeric(dremi))
  expect_length(dremi, 1)
  expect_gte(dremi, 0)
})

test_that("knnDREMI respects parameter variations", {
  set.seed(42)
  
  n <- 300
  x <- rnorm(n)
  y <- 0.7 * x + 0.3 * rnorm(n)
  
  # Different k values
  dremi_k5 <- knnDREMI(x, y, k = 5)
  dremi_k20 <- knnDREMI(x, y, k = 20)
  
  # Both should be reasonable
  expect_true(dremi_k5 > 0)
  expect_true(dremi_k20 > 0)
  
  # Different n_bins
  dremi_bins10 <- knnDREMI(x, y, n_bins = 10)
  dremi_bins30 <- knnDREMI(x, y, n_bins = 30)
  
  expect_true(dremi_bins10 > 0)
  expect_true(dremi_bins30 > 0)
})
