# Tests for knnDREMI and related functions

test_that("knnDREMI computes correctly", {
  set.seed(42)
  
  # Create correlated data
  n <- 200
  x <- rnorm(n)
  y <- 0.8 * x + 0.2 * rnorm(n)  # Highly correlated
  
  dremi <- knnDREMI(x, y, k = 10)
  
  # DREMI should be between 0 and 1

  expect_gte(dremi, 0)
  expect_lte(dremi, 1)
  
  # For highly correlated data, DREMI should be relatively high
  expect_gt(dremi, 0.3)
})

test_that("knnDREMI handles uncorrelated data", {
  set.seed(42)
  
  # Uncorrelated data
  n <- 200
  x <- rnorm(n)
  y <- rnorm(n)  # Independent
  
  dremi <- knnDREMI(x, y, k = 10)
  
  # DREMI should be low for independent data
  expect_gte(dremi, 0)
  expect_lte(dremi, 1)
  expect_lt(dremi, 0.3)  # Should be relatively low
})

test_that("knnDREMI returns DREVI when requested", {
  set.seed(42)
  
  n <- 200
  x <- rnorm(n)
  y <- 0.8 * x + 0.2 * rnorm(n)
  
  result <- knnDREMI(x, y, k = 10, return_drevi = TRUE)
  
  expect_type(result, "list")
  expect_true("dremi" %in% names(result))
  expect_true("drevi" %in% names(result))
  
  expect_gte(result$dremi, 0)
  expect_lte(result$dremi, 1)
  expect_gte(result$drevi, 0)
  expect_lte(result$drevi, 1)
})

test_that("knnDREMI validates inputs", {
  x <- rnorm(100)
  y <- rnorm(50)  # Different length
  
  expect_error(knnDREMI(x, y), "same length")
})

test_that("magic_knnDREMI works with magic result", {
  skip_if_not_installed("MAGICR")
  
  data(magic_testdata)
  
  result <- magic(magic_testdata, t = 2, verbose = FALSE)
  
  # Get two genes
  gene1 <- colnames(magic_testdata)[1]
  gene2 <- colnames(magic_testdata)[2]
  
  dremi <- magic_knnDREMI(result, gene1, gene2)
  
  expect_gte(dremi, 0)
  expect_lte(dremi, 1)
})

test_that("knnDREMI handles edge cases", {
  # Very small k
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)
  
  dremi <- knnDREMI(x, y, k = 3)
  expect_gte(dremi, 0)
  expect_lte(dremi, 1)
  
  # Large k
  dremi <- knnDREMI(x, y, k = 50)
  expect_gte(dremi, 0)
  expect_lte(dremi, 1)
})

test_that("knnDREMI is asymmetric", {
  set.seed(42)
  
  # Create asymmetric relationship
  n <- 200
  x <- rnorm(n)
  y <- x^2 + 0.1 * rnorm(n)  # y depends on x, but not vice versa
  
  dremi_xy <- knnDREMI(x, y, k = 10)
  dremi_yx <- knnDREMI(y, x, k = 10)
  
  # Should give different values
  expect_false(abs(dremi_xy - dremi_yx) < 0.01)
})
