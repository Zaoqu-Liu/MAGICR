#' @title MAGIC: Markov Affinity-based Graph Imputation of Cells
#' @description Main function for running MAGIC on single-cell data.
#' @name magic
NULL

#' Run MAGIC Imputation
#'
#' @description Performs MAGIC (Markov Affinity-based Graph Imputation of Cells)
#'   for denoising and imputation of single-cell data.
#'
#' @param data Input data. Can be a matrix, sparse Matrix (dgCMatrix), data.frame,
#'   Seurat object, or SingleCellExperiment object.
#' @param genes Character vector of gene names, integer vector of gene indices,
#'   or one of the special values "all_genes" or "pca_only". If NULL, returns
#'   all genes. Default is NULL.
#' @param knn Integer. Number of nearest neighbors for kernel bandwidth.
#'   Default is 5.
#' @param knn_max Integer or NULL. Maximum number of neighbors to consider.
#'   If NULL, defaults to 3 * knn. Default is NULL.
#' @param decay Numeric or NULL. Alpha parameter for kernel decay rate.
#'   Higher values create sharper kernels. If NULL, uses unweighted kNN.
#'   Default is 1.
#' @param t Integer or "auto". Number of diffusion steps. If "auto", optimal t
#'   is selected automatically using Procrustes analysis. Default is 3.
#' @param npca Integer or NULL. Number of principal components for distance
#'   calculation. If NULL or greater than number of genes, PCA is not performed.
#'   Default is 100.
#' @param solver Character. Either "exact" (compute on all genes) or
#'   "approximate" (compute in PCA space, then project back). The "approximate"
#'   solver is faster but may return negative values. Default is "exact".
#' @param init A previous MAGIC result object. If provided, reuses the
#'   diffusion operator from that run. Default is NULL.
#' @param t_max Integer. Maximum t to test when t = "auto". Default is 20.
#' @param distance Character. Distance metric for kNN. One of "euclidean",
#'   "cosine", "manhattan", or "correlation". Default is "euclidean".
#' @param n_jobs Integer. Number of parallel workers. Use -1 for all available.
#'   Default is 1 (no parallelization).
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#' @param verbose Logical or integer. Verbosity level. TRUE or 1 for basic
#'   messages, 2 for detailed messages, 0 or FALSE for silent. Default is TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return For matrices/data.frames: A "magic" object containing:
#'   \item{result}{Matrix of imputed expression values}
#'   \item{diff_op}{The diffusion operator used}
#'   \item{params}{List of parameters used}
#'   \item{pca}{PCA results if npca was used, NULL otherwise}
#'
#'   For Seurat objects: The input object with a new assay containing imputed values.
#'
#'   For SingleCellExperiment: The input object with imputed values in a new assay.
#'
#' @details MAGIC uses diffusion on a cell-cell similarity graph to denoise
#'   single-cell RNA-seq data. The algorithm:
#'   \enumerate{
#'     \item Optionally reduces dimensions via PCA
#'     \item Builds a k-nearest neighbor graph
#'     \item Constructs an alpha-decaying kernel
#'     \item Creates a row-stochastic diffusion operator

#'     \item Applies the diffusion operator t times to smooth expression
#'   }
#'
#'   The level of smoothing is controlled by the parameter t. Higher values
#'   of t result in more smoothing but may over-smooth the data.
#'
#' @references
#' van Dijk, D., Sharma, R., Nainys, J., Yim, K., Kathail, P., Carr, A.J.,
#' Burdziak, C., Moon, K.R., Chaffer, C.L., Pattabiraman, D., Bierie, B.,
#' Mazutis, L., Wolf, G., Krishnaswamy, S., & Pe'er, D. (2018).
#' Recovering gene interactions from single-cell data using data diffusion.
#' Cell, 174(3), 716-729.e27. \doi{10.1016/j.cell.2018.05.061}
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic usage with matrix
#' set.seed(42)
#' data <- matrix(rpois(5000, lambda = 2), nrow = 500, ncol = 10)
#' colnames(data) <- paste0("Gene", 1:10)
#' rownames(data) <- paste0("Cell", 1:500)
#'
#' # Run MAGIC
#' result <- magic(data, t = 3, verbose = TRUE)
#'
#' # Access imputed data
#' imputed <- as.matrix(result)
#'
#' # Run on specific genes
#' result_subset <- magic(data, genes = c("Gene1", "Gene2", "Gene3"))
#'
#' # With automatic t selection
#' result_auto <- magic(data, t = "auto")
#'
#' # With Seurat object
#' library(Seurat)
#' seurat_obj <- CreateSeuratObject(counts = t(data))
#' seurat_obj <- NormalizeData(seurat_obj)
#' seurat_obj <- magic(seurat_obj)
#' }
magic <- function(data, ...) {
  UseMethod("magic")
}

#' @rdname magic
#' @export
magic.default <- function(data,
                          genes = NULL,
                          knn = 5,
                          knn_max = NULL,
                          decay = 1,
                          t = 3,
                          npca = 100,
                          solver = c("exact", "approximate"),
                          init = NULL,
                          t_max = 20,
                          distance = c("euclidean", "cosine", "manhattan", "correlation"),
                          n_jobs = 1,
                          seed = NULL,
                          verbose = TRUE,
                          ...) {
  
  # Convert to matrix if needed
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  magic.matrix(
    data = data,
    genes = genes,
    knn = knn,
    knn_max = knn_max,
    decay = decay,
    t = t,
    npca = npca,
    solver = solver,
    init = init,
    t_max = t_max,
    distance = distance,
    n_jobs = n_jobs,
    seed = seed,
    verbose = verbose,
    ...
  )
}

#' @rdname magic
#' @export
magic.matrix <- function(data,
                         genes = NULL,
                         knn = 5,
                         knn_max = NULL,
                         decay = 1,
                         t = 3,
                         npca = 100,
                         solver = c("exact", "approximate"),
                         init = NULL,
                         t_max = 20,
                         distance = c("euclidean", "cosine", "manhattan", "correlation"),
                         n_jobs = 1,
                         seed = NULL,
                         verbose = TRUE,
                         ...) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Validate parameters
  solver <- match.arg(solver)
  distance <- match.arg(distance)
  
  check_positive_int(knn, "knn")
  if (!is.null(knn_max)) check_positive_int(knn_max, "knn_max")
  if (!is.null(npca)) check_positive_int(npca, "npca")
  check_positive_int(t_max, "t_max")
  if (!identical(t, "auto")) check_positive_int(t, "t")
  if (!is.null(decay)) check_positive(decay, "decay")
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
    vmessage(sprintf("Random seed set to %d", seed), verbose, level = 2)
  }
  
  # Setup parallel backend if requested
  if (n_jobs != 1) {
    setup_parallel(n_jobs, verbose)
  }
  
  # Validate and prepare data
  vmessage(sprintf("Running MAGIC on %d cells and %d genes...",
                   nrow(data), ncol(data)), verbose)
  
  data <- validate_data(data, verbose)
  
  # Store original data info
  original_rownames <- safe_rownames(data)
  original_colnames <- safe_colnames(data)
  
  # Parse genes argument
  genes_parsed <- parse_genes(genes, original_colnames)
  return_pca_only <- identical(genes_parsed$mode, "pca_only")
  gene_indices <- genes_parsed$indices
  
  # Check if we can reuse previous diffusion operator
  diff_op <- NULL
  pca_result <- NULL
  
  if (!is.null(init) && inherits(init, "magic")) {
    vmessage("  Using diffusion operator from previous run...", verbose)
    diff_op <- init$diff_op
    pca_result <- init$pca
    
    # Validate dimensions match
    if (!is.null(diff_op) && nrow(diff_op) != nrow(data)) {
      warning("Dimensions don't match previous run. Recomputing diffusion operator.")
      diff_op <- NULL
      pca_result <- NULL
    }
  }
  
  # Compute diffusion operator if needed
  if (is.null(diff_op)) {
    # PCA reduction if requested
    if (!is.null(npca) && ncol(data) > npca) {
      vmessage(sprintf("  Performing PCA (n_components=%d)...", npca), verbose)
      pca_result <- timed(
        compute_pca(data, npca, seed),
        "PCA", verbose
      )
      data_for_graph <- pca_result$x
    } else {
      data_for_graph <- data
      vmessage("  Skipping PCA (fewer genes than npca)", verbose, level = 2)
    }
    
    # Build diffusion operator
    vmessage("  Building diffusion operator...", verbose)
    diff_op <- timed(
      build_diffusion_operator(
        data = data_for_graph,
        knn = knn,
        knn_max = knn_max,
        decay = decay,
        distance = distance,
        verbose = verbose
      ),
      "Diffusion operator", verbose
    )
  }
  
  # Perform imputation
  if (solver == "approximate" || return_pca_only) {
    # Approximate solver: impute in PCA space
    if (is.null(pca_result)) {
      vmessage("  Computing PCA for approximate solver...", verbose)
      pca_result <- compute_pca(data, npca %||% 100, seed)
    }
    
    vmessage("  Imputing in PCA space...", verbose)
    imputed_pca <- timed(
      magic_impute(
        data = pca_result$x,
        diff_op = diff_op,
        t = t,
        t_max = t_max,
        threshold = 0.001,
        solver = solver,
        verbose = verbose
      ),
      "Imputation", verbose
    )
    
    if (return_pca_only) {
      imputed <- imputed_pca
      output_colnames <- paste0("PC", seq_len(ncol(imputed_pca)))
    } else {
      # Project back to gene space
      vmessage("  Projecting back to gene space...", verbose)
      imputed <- project_from_pca(imputed_pca, pca_result, gene_indices)
      output_colnames <- if (is.null(gene_indices)) {
        original_colnames
      } else {
        original_colnames[gene_indices]
      }
    }
  } else {
    # Exact solver: impute in gene space
    if (!is.null(gene_indices)) {
      data_to_impute <- data[, gene_indices, drop = FALSE]
    } else {
      data_to_impute <- data
    }
    
    vmessage("  Imputing gene expression...", verbose)
    imputed <- timed(
      magic_impute(
        data = data_to_impute,
        diff_op = diff_op,
        t = t,
        t_max = t_max,
        threshold = 0.001,
        solver = solver,
        verbose = verbose
      ),
      "Imputation", verbose
    )
    
    output_colnames <- if (is.null(gene_indices)) {
      original_colnames
    } else {
      original_colnames[gene_indices]
    }
  }
  
  # Set names
  rownames(imputed) <- original_rownames
  colnames(imputed) <- output_colnames
  
  # Get optimal t if auto was used
  optimal_t <- attr(imputed, "optimal_t")
  if (is.null(optimal_t)) optimal_t <- t
  
  # Create result object
  result <- list(
    result = as.data.frame(imputed),
    diff_op = diff_op,
    params = list(
      genes = genes,
      knn = knn,
      knn_max = knn_max %||% (knn * 3),
      decay = decay,
      t = optimal_t,
      npca = npca,
      solver = solver,
      distance = distance,
      n_cells = nrow(data),
      n_genes = ncol(data)
    ),
    pca = pca_result
  )
  
  class(result) <- c("magic", "list")
  
  # Report timing
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  vmessage(sprintf("MAGIC complete in %s", format_time(elapsed)), verbose)
  
  return(result)
}

#' @rdname magic
#' @export
magic.dgCMatrix <- function(data, ...) {
  magic.matrix(as.matrix(data), ...)
}

#' @rdname magic
#' @export
magic.data.frame <- function(data, ...) {
  magic.matrix(as.matrix(data), ...)
}

#' Parse Genes Argument
#' @param genes User-provided genes argument
#' @param colnames Column names of data
#' @return List with mode and indices
#' @keywords internal
parse_genes <- function(genes, colnames) {
  if (is.null(genes) || (length(genes) == 1 && genes == "all_genes")) {
    return(list(mode = "all", indices = NULL))
  }
  
  if (length(genes) == 1 && genes == "pca_only") {
    return(list(mode = "pca_only", indices = NULL))
  }
  
  if (is.character(genes)) {
    # Gene names
    missing <- genes[!genes %in% colnames]
    if (length(missing) > 0) {
      warning(sprintf("Genes not found: %s", paste(missing, collapse = ", ")))
    }
    indices <- which(colnames %in% genes)
    if (length(indices) == 0) {
      stop("No valid genes found", call. = FALSE)
    }
    return(list(mode = "subset", indices = indices))
  }
  
  if (is.numeric(genes)) {
    # Gene indices
    genes <- as.integer(genes)
    if (any(genes < 1) || any(genes > length(colnames))) {
      stop("Gene indices out of range", call. = FALSE)
    }
    return(list(mode = "subset", indices = genes))
  }
  
  stop("Invalid 'genes' argument", call. = FALSE)
}

#' Compute PCA
#' @param data Data matrix
#' @param npca Number of components
#' @param seed Random seed
#' @return PCA result with x (scores) and rotation (loadings)
#' @keywords internal
compute_pca <- function(data, npca, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Use irlba for fast PCA
  npca <- min(npca, min(dim(data)) - 1)
  
  # Center data
  data_centered <- scale(data, center = TRUE, scale = FALSE)
  
  # Fast randomized SVD via irlba
  pca_result <- irlba::prcomp_irlba(data_centered, n = npca, center = FALSE)
  
  return(list(
    x = pca_result$x,
    rotation = pca_result$rotation,
    center = attr(data_centered, "scaled:center"),
    sdev = pca_result$sdev
  ))
}

#' Project from PCA Space Back to Gene Space
#' @param imputed_pca Imputed PCA scores
#' @param pca_result PCA result object
#' @param gene_indices Indices of genes to return (NULL for all)
#' @return Matrix in gene space
#' @keywords internal
project_from_pca <- function(imputed_pca, pca_result, gene_indices = NULL) {
  # X_reconstructed = X_pca %*% t(rotation) + center
  result <- imputed_pca %*% t(pca_result$rotation)
  result <- sweep(result, 2, pca_result$center, "+")
  
  if (!is.null(gene_indices)) {
    result <- result[, gene_indices, drop = FALSE]
  }
  
  return(result)
}

#' Setup Parallel Backend
#' @param n_jobs Number of workers
#' @param verbose Verbosity level
#' @keywords internal
setup_parallel <- function(n_jobs, verbose = FALSE) {
  if (n_jobs == -1) {
    n_jobs <- future::availableCores()
  }
  
  if (n_jobs > 1) {
    vmessage(sprintf("  Setting up parallel backend with %d workers...", n_jobs),
             verbose)
    future::plan(future::multisession, workers = n_jobs)
  }
  
  invisible(NULL)
}

# S3 Methods --------------------------------------------------------------

#' @rdname magic
#' @export
print.magic <- function(x, ...) {
  cat("MAGIC Result\n")
  cat("============\n")
  cat(sprintf("  Cells: %d\n", nrow(x$result)))
  cat(sprintf("  Genes: %d\n", ncol(x$result)))
  cat(sprintf("  Parameters:\n"))
  cat(sprintf("    knn: %d\n", x$params$knn))
  cat(sprintf("    decay: %s\n", x$params$decay %||% "NULL"))
  cat(sprintf("    t: %d\n", x$params$t))
  cat(sprintf("    solver: %s\n", x$params$solver))
  cat(sprintf("    distance: %s\n", x$params$distance))
  if (!is.null(x$params$npca)) {
    cat(sprintf("    npca: %d\n", x$params$npca))
  }
  invisible(x)
}

#' @rdname magic
#' @export
summary.magic <- function(object, ...) {
  summary(object$result)
}

#' @rdname magic
#' @export
as.matrix.magic <- function(x, ...) {
  as.matrix(x$result)
}

#' @rdname magic
#' @export
as.data.frame.magic <- function(x, ...) {
  x$result
}

#' Subset MAGIC Result
#' @param x MAGIC object
#' @param i Row indices
#' @param j Column indices
#' @param ... Additional arguments
#' @export
`[.magic` <- function(x, i, j, ...) {
  if (missing(i) && missing(j)) {
    return(x$result)
  }
  if (missing(i)) {
    return(x$result[, j, drop = FALSE])
  }
  if (missing(j)) {
    return(x$result[i, , drop = FALSE])
  }
  x$result[i, j, drop = FALSE]
}
