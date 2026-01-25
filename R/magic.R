#' @title MAGIC: Markov Affinity-based Graph Imputation of Cells
#' @description Main function for MAGIC imputation of single-cell data.
#' @name magic
NULL

#' Run MAGIC
#'
#' @description Performs MAGIC (Markov Affinity-based Graph Imputation of Cells)
#'   for denoising and imputation of single-cell data.
#'
#' @param data Input data. Matrix, sparse Matrix, data.frame, Seurat, or
#'   SingleCellExperiment object.
#' @param genes Gene names, indices, "all_genes", or "pca_only". Default is NULL
#'   (all genes).
#' @param knn Integer. Neighbors for kernel bandwidth. Default is 5.
#' @param knn_max Integer or NULL. Maximum neighbors. Default is 3 * knn.
#' @param decay Numeric or NULL. Alpha decay parameter. Default is 1.
#' @param t Integer or "auto". Diffusion steps. Default is 3.
#' @param npca Integer or NULL. PCA components. Default is 100.
#' @param solver "exact" or "approximate". Default is "exact".
#' @param init Previous MAGIC result. Default is NULL.
#' @param t_max Integer. Maximum t for automatic selection. Default is 20.
#' @param knn_dist Distance metric. Default is "euclidean".
#' @param n_jobs Integer. Parallel workers (-1 for all). Default is 1.
#' @param seed Integer or NULL. Random seed. Default is NULL.
#' @param verbose Logical or integer. Default is TRUE.
#' @param ... Additional arguments.
#'
#' @return A "magic" object (for matrices) or modified input object (for Seurat/SCE).
#'
#' @details MAGIC denoises single-cell data using diffusion on a cell similarity graph:
#'   \enumerate{
#'     \item Optional PCA dimension reduction
#'     \item Build k-nearest neighbor graph
#'     \item Construct alpha-decaying kernel
#'     \item Create row-stochastic diffusion operator
#'     \item Apply diffusion operator t times
#'   }
#'   
#'   The imputation formula is: \eqn{X_{imputed} = P^t X}
#'
#' @references
#' van Dijk, D., et al. (2018). Recovering gene interactions from single-cell
#' data using data diffusion. Cell, 174(3), 716-729.
#' \doi{10.1016/j.cell.2018.05.061}
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic usage
#' set.seed(42)
#' data <- matrix(rpois(5000, 2), nrow = 500, ncol = 10)
#' result <- magic(data, t = 3)
#' imputed <- as.matrix(result)
#'
#' # Automatic t selection
#' result <- magic(data, t = "auto")
#'
#' # Specific genes
#' result <- magic(data, genes = c(1, 2, 3))
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
                          knn_dist = c("euclidean", "cosine", "manhattan", "correlation"),
                          n_jobs = 1,
                          seed = NULL,
                          verbose = TRUE,
                          ...) {
  
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
    knn_dist = knn_dist,
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
                         knn_dist = c("euclidean", "cosine", "manhattan", "correlation"),
                         n_jobs = 1,
                         seed = NULL,
                         verbose = TRUE,
                         ...) {
  
  start_time <- Sys.time()
  
  # Validate parameters
  solver <- match.arg(solver)
  knn_dist <- match.arg(knn_dist)
  
  check_positive_int(knn, "knn")
  if (!is.null(knn_max)) check_positive_int(knn_max, "knn_max")
  if (!is.null(npca)) check_positive_int(npca, "npca")
  check_positive_int(t_max, "t_max")
  if (!identical(t, "auto")) check_positive_int(t, "t")
  if (!is.null(decay)) check_positive(decay, "decay")
  
  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Parallel setup
  if (n_jobs != 1) {
    setup_parallel(n_jobs, verbose)
  }
  
  # Validate data
  vmessage(sprintf("Running MAGIC on %d cells and %d genes...",
                   nrow(data), ncol(data)), verbose)
  
  data <- validate_data(data, verbose)
  
  # Store names
  original_rownames <- safe_rownames(data)
  original_colnames <- safe_colnames(data)
  
  # Parse genes
  genes_parsed <- parse_genes(genes, original_colnames)
  return_pca_only <- identical(genes_parsed$mode, "pca_only")
  gene_indices <- genes_parsed$indices
  
  # Reuse diffusion operator if provided
  diff_op <- NULL
  pca_result <- NULL
  
  if (!is.null(init) && inherits(init, "magic")) {
    vmessage("  Using previous diffusion operator...", verbose)
    diff_op <- init$diff_op
    pca_result <- init$pca
    
    if (!is.null(diff_op) && nrow(diff_op) != nrow(data)) {
      warning("Dimension mismatch. Recomputing diffusion operator.")
      diff_op <- NULL
      pca_result <- NULL
    }
  }
  
  # Build diffusion operator
  if (is.null(diff_op)) {
    # PCA
    if (!is.null(npca) && ncol(data) > npca) {
      vmessage(sprintf("  PCA (n=%d)...", npca), verbose)
      pca_result <- timed(compute_pca(data, npca, seed), "PCA", verbose)
      data_for_graph <- pca_result$x
    } else {
      data_for_graph <- data
    }
    
    # Diffusion operator
    vmessage("  Building diffusion operator...", verbose)
    diff_op <- timed(
      build_diffusion_operator(
        data = data_for_graph,
        knn = knn,
        knn_max = knn_max,
        decay = decay,
        distance = knn_dist,
        verbose = verbose
      ),
      "Diffusion operator", verbose
    )
  }
  
  # Imputation
  if (solver == "approximate" || return_pca_only) {
    if (is.null(pca_result)) {
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
      vmessage("  Projecting to gene space...", verbose)
      imputed <- project_from_pca(imputed_pca, pca_result, gene_indices)
      output_colnames <- if (is.null(gene_indices)) {
        original_colnames
      } else {
        original_colnames[gene_indices]
      }
    }
  } else {
    # Exact solver
    if (!is.null(gene_indices)) {
      data_to_impute <- data[, gene_indices, drop = FALSE]
    } else {
      data_to_impute <- data
    }
    
    vmessage("  Imputing...", verbose)
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
  
  # Get t used
  optimal_t <- attr(imputed, "optimal_t")
  if (is.null(optimal_t)) optimal_t <- t
  
  # Create result
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
      knn_dist = knn_dist,
      n_cells = nrow(data),
      n_genes = ncol(data)
    ),
    pca = pca_result
  )
  
  class(result) <- c("magic", "list")
  
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

#' Parse genes argument
#' @keywords internal
parse_genes <- function(genes, colnames) {
  if (is.null(genes) || (length(genes) == 1 && genes == "all_genes")) {
    return(list(mode = "all", indices = NULL))
  }
  
  if (length(genes) == 1 && genes == "pca_only") {
    return(list(mode = "pca_only", indices = NULL))
  }
  
  if (is.character(genes)) {
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
    genes <- as.integer(genes)
    if (any(genes < 1) || any(genes > length(colnames))) {
      stop("Gene indices out of range", call. = FALSE)
    }
    return(list(mode = "subset", indices = genes))
  }
  
  stop("Invalid 'genes' argument", call. = FALSE)
}

#' Compute PCA
#' @keywords internal
compute_pca <- function(data, npca, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  npca <- min(npca, min(dim(data)) - 1)
  data_centered <- scale(data, center = TRUE, scale = FALSE)
  pca_result <- irlba::prcomp_irlba(data_centered, n = npca, center = FALSE)
  
  return(list(
    x = pca_result$x,
    rotation = pca_result$rotation,
    center = attr(data_centered, "scaled:center"),
    sdev = pca_result$sdev
  ))
}

#' Project from PCA to gene space
#' @keywords internal
project_from_pca <- function(imputed_pca, pca_result, gene_indices = NULL) {
  result <- imputed_pca %*% t(pca_result$rotation)
  result <- sweep(result, 2, pca_result$center, "+")
  
  if (!is.null(gene_indices)) {
    result <- result[, gene_indices, drop = FALSE]
  }
  
  return(result)
}

#' Setup parallel backend
#' @keywords internal
setup_parallel <- function(n_jobs, verbose = FALSE) {
  if (n_jobs == -1) {
    n_jobs <- future::availableCores()
  }
  
  if (n_jobs > 1) {
    vmessage(sprintf("  Parallel: %d workers", n_jobs), verbose)
    future::plan(future::multisession, workers = n_jobs)
  }
  
  invisible(NULL)
}

# S3 Methods ----------------------------------------------------------------

#' @method print magic
#' @export
print.magic <- function(x, ...) {
  cat("MAGIC Result\n")
  cat("============\n")
  cat(sprintf("  Cells: %d\n", nrow(x$result)))
  cat(sprintf("  Genes: %d\n", ncol(x$result)))
  cat("  Parameters:\n")
  cat(sprintf("    knn: %d\n", x$params$knn))
  cat(sprintf("    decay: %s\n", x$params$decay %||% "NULL"))
  cat(sprintf("    t: %s\n", x$params$t))
  cat(sprintf("    solver: %s\n", x$params$solver))
  cat(sprintf("    knn_dist: %s\n", x$params$knn_dist))
  if (!is.null(x$params$npca)) {
    cat(sprintf("    npca: %d\n", x$params$npca))
  }
  invisible(x)
}

#' @method summary magic
#' @export
summary.magic <- function(object, ...) {
  summary(object$result)
}

#' @method as.matrix magic
#' @export
as.matrix.magic <- function(x, ...) {
  as.matrix(x$result)
}

#' @method as.data.frame magic
#' @export
as.data.frame.magic <- function(x, ...) {
  x$result
}

#' Subset MAGIC result
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
