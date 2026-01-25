// MAGICR - C++ acceleration functions
// Author: Zaoqu Liu
// This file contains optimized C++ implementations of performance-critical
// operations for the MAGICR package.
// Note: No OpenMP used to avoid conflicts on different systems.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Check if Rcpp functions are available
//' @keywords internal
// [[Rcpp::export]]
bool rcpp_available() {
  return true;
}

//' Fast Euclidean Distance Computation
//' 
//' Computes pairwise Euclidean distances efficiently using matrix operations.
//' @param X Matrix where rows are observations
//' @return Distance matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_euclidean_dist(const arma::mat& X) {
  int n = X.n_rows;
  arma::mat D(n, n);
  
  // Compute squared norms
  arma::vec sq_norms = sum(X % X, 1);
  
  // D_ij^2 = ||x_i||^2 + ||x_j||^2 - 2 * x_i . x_j
  arma::mat XXT = X * X.t();
  
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      double d2 = sq_norms(i) + sq_norms(j) - 2.0 * XXT(i, j);
      d2 = std::max(0.0, d2);  // Numerical stability
      D(i, j) = std::sqrt(d2);
      D(j, i) = D(i, j);
    }
  }
  
  return D;
}

//' Fast Alpha-Decaying Kernel Computation
//' 
//' Computes the alpha-decaying kernel values from distances and bandwidths.
//' @param distances Matrix of distances (n x k)
//' @param epsilon Vector of adaptive bandwidths (n)
//' @param alpha Decay parameter
//' @param threshold Threshold for zeroing small values
//' @return Kernel values matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_alpha_kernel(const arma::mat& distances,
                           const arma::vec& epsilon,
                           double alpha,
                           double threshold) {
  int n = distances.n_rows;
  int k = distances.n_cols;
  
  arma::mat K(n, k);
  
  for (int i = 0; i < n; i++) {
    double eps_i = epsilon(i);
    if (eps_i <= 0) eps_i = 1e-10;  // Avoid division by zero
    
    for (int j = 0; j < k; j++) {
      double d_normalized = distances(i, j) / eps_i;
      double kernel_val = std::exp(-std::pow(d_normalized, alpha));
      
      if (kernel_val < threshold) {
        kernel_val = 0.0;
      }
      K(i, j) = kernel_val;
    }
  }
  
  return K;
}

//' Fast Sparse Matrix-Dense Matrix Multiplication
//' 
//' Optimized multiplication of sparse diffusion operator with dense data matrix.
//' This is the core operation in MAGIC imputation.
//' @param P_i Row indices of sparse matrix (0-indexed)
//' @param P_p Column pointers of sparse matrix (0-indexed)
//' @param P_x Values of sparse matrix
//' @param X Dense matrix to multiply
//' @param nrow Number of rows in P
//' @return Result of P %*% X
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_sparse_dense_mult(const arma::uvec& P_i,
                                 const arma::uvec& P_p,
                                 const arma::vec& P_x,
                                 const arma::mat& X,
                                 int nrow) {
  int ncol_X = X.n_cols;
  arma::mat result(nrow, ncol_X, arma::fill::zeros);
  
  // CSC format: P_p[j] to P_p[j+1]-1 are the row indices for column j
  int ncol_P = P_p.n_elem - 1;
  
  for (int j = 0; j < ncol_P; j++) {
    for (unsigned int idx = P_p(j); idx < P_p(j + 1); idx++) {
      int i = P_i(idx);
      double val = P_x(idx);
      // result[i, :] += val * X[j, :]
      result.row(i) += val * X.row(j);
    }
  }
  
  return result;
}

//' Fast Matrix Power
//' 
//' Computes matrix power using binary exponentiation.
//' @param M Square matrix
//' @param n Power (positive integer)
//' @return M^n
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_matrix_power(const arma::mat& M, int n) {
  int size = M.n_rows;
  
  if (n == 0) {
    return arma::eye(size, size);
  }
  if (n == 1) {
    return M;
  }
  
  arma::mat result = arma::eye(size, size);
  arma::mat base = M;
  
  while (n > 0) {
    if (n % 2 == 1) {
      result = result * base;
    }
    base = base * base;
    n = n / 2;
  }
  
  return result;
}

//' Fast Procrustes Disparity
//' 
//' Computes Procrustes disparity between two matrices.
//' @param X First matrix (n x p)
//' @param Y Second matrix (n x p)
//' @return Disparity value
//' @keywords internal
// [[Rcpp::export]]
double cpp_procrustes_disparity(const arma::mat& X, const arma::mat& Y) {
  // Center matrices
  arma::mat X_centered = X - arma::repmat(arma::mean(X, 0), X.n_rows, 1);
  arma::mat Y_centered = Y - arma::repmat(arma::mean(Y, 0), Y.n_rows, 1);
  
  // Frobenius norms
  double X_norm = arma::norm(X_centered, "fro");
  double Y_norm = arma::norm(Y_centered, "fro");
  
  if (X_norm < 1e-10 || Y_norm < 1e-10) {
    return 0.0;
  }
  
  // Normalize
  X_centered /= X_norm;
  Y_centered /= Y_norm;
  
  // SVD of X'Y
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, X_centered.t() * Y_centered);
  
  // Optimal rotation
  arma::mat R = V * U.t();
  
  // Disparity
  arma::mat Y_aligned = Y_centered * R;
  double disparity = arma::accu(arma::square(X_centered - Y_aligned));
  
  return disparity;
}

//' Row Normalize Sparse Matrix
//' 
//' Makes a sparse matrix row-stochastic (rows sum to 1).
//' @param i Row indices (1-indexed, R style)
//' @param j Column indices (1-indexed, R style)
//' @param x Values
//' @param nrow Number of rows
//' @param ncol Number of columns
//' @return List with normalized values and row sums
//' @keywords internal
// [[Rcpp::export]]
List cpp_row_normalize(const arma::uvec& i,
                       const arma::uvec& j,
                       const arma::vec& x,
                       int nrow,
                       int ncol) {
  
  // Compute row sums
  arma::vec row_sums(nrow, arma::fill::zeros);
  for (unsigned int idx = 0; idx < i.n_elem; idx++) {
    row_sums(i(idx) - 1) += x(idx);  // R is 1-indexed
  }
  
  // Normalize values
  arma::vec x_norm(x.n_elem);
  for (unsigned int idx = 0; idx < i.n_elem; idx++) {
    int row = i(idx) - 1;
    if (row_sums(row) > 0) {
      x_norm(idx) = x(idx) / row_sums(row);
    } else {
      x_norm(idx) = 0.0;
    }
  }
  
  return List::create(
    Named("x") = x_norm,
    Named("row_sums") = row_sums
  );
}

//' Fast KNN Bandwidth Computation
//' 
//' Extracts the k-th nearest neighbor distance for adaptive bandwidth.
//' @param distances Matrix of distances (n x knn_max)
//' @param k Which neighbor to use (1-indexed)
//' @return Vector of bandwidths
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_get_bandwidth(const arma::mat& distances, int k) {
  int n = distances.n_rows;
  arma::vec epsilon(n);
  
  // k is 1-indexed, so column index is k-1
  // But note that column 0 might be the point itself (distance 0)
  // So we use column k (the k-th neighbor after self)
  int col_idx = std::min(k, (int)distances.n_cols - 1);
  
  for (int i = 0; i < n; i++) {
    epsilon(i) = distances(i, col_idx);
    if (epsilon(i) <= 0) {
      epsilon(i) = 1e-10;  // Avoid zero bandwidth
    }
  }
  
  return epsilon;
}

//' Dense Matrix Imputation Step
//' 
//' Performs one step of diffusion: X_new = P %*% X
//' Optimized for dense diffusion operator.
//' @param P Diffusion operator (dense)
//' @param X Data matrix
//' @return Diffused data
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_diffuse_step(const arma::mat& P, const arma::mat& X) {
  return P * X;
}

//' Multiple Diffusion Steps
//' 
//' Performs t steps of diffusion efficiently.
//' @param P Diffusion operator
//' @param X Data matrix
//' @param t Number of steps
//' @return Diffused data after t steps
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_diffuse_t_steps(const arma::mat& P, const arma::mat& X, int t) {
  if (t <= 0) return X;
  
  arma::mat result = X;
  for (int i = 0; i < t; i++) {
    result = P * result;
  }
  
  return result;
}

//' Automatic t Selection with Procrustes
//' 
//' Iteratively applies diffusion and computes Procrustes disparity.
//' @param P Diffusion operator
//' @param X Data matrix
//' @param t_max Maximum t to try
//' @param threshold Convergence threshold
//' @param subsample_cols Column indices to subsample (0-indexed), empty for all
//' @return List with imputed data, optimal t, and error vector
//' @keywords internal
// [[Rcpp::export]]
List cpp_auto_t_impute(const arma::mat& P,
                       const arma::mat& X,
                       int t_max,
                       double threshold,
                       const arma::uvec& subsample_cols) {
  
  arma::mat X_curr = X;
  arma::mat X_prev;
  arma::vec error_vec(t_max);
  int t_opt = t_max;
  
  // Prepare subsampled view
  arma::mat X_sub_curr, X_sub_prev;
  bool use_subsample = subsample_cols.n_elem > 0;
  
  for (int t = 0; t < t_max; t++) {
    X_prev = X_curr;
    X_curr = P * X_curr;
    
    // Compute error on subsampled columns
    if (use_subsample) {
      X_sub_prev = X_prev.cols(subsample_cols);
      X_sub_curr = X_curr.cols(subsample_cols);
      error_vec(t) = cpp_procrustes_disparity(X_sub_prev, X_sub_curr);
    } else {
      error_vec(t) = cpp_procrustes_disparity(X_prev, X_curr);
    }
    
    // Check convergence
    if (error_vec(t) < threshold && t_opt == t_max) {
      t_opt = t + 1;  // 1-indexed
      break;
    }
  }
  
  return List::create(
    Named("imputed") = X_curr,
    Named("t_opt") = t_opt,
    Named("error_vec") = error_vec.head(t_opt)
  );
}
