#' @title Seurat Integration
#' @description Functions for MAGIC integration with Seurat objects (v4 and v5).
#' @name seurat
NULL

#' Detect Seurat Version
#' @return Integer version number (4 or 5)
#' @keywords internal
detect_seurat_version <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required but not installed.", call. = FALSE)
  }
  
  version <- utils::packageVersion("Seurat")
  major <- as.integer(version[[1]][1])
  
  return(major)
}

#' Get Data from Seurat Object
#' @description Extracts normalized data from Seurat object, compatible with v4 and v5.
#' @param object Seurat object
#' @param assay Assay name
#' @param layer Layer/slot name
#' @return Matrix (genes x cells)
#' @keywords internal
get_seurat_data <- function(object, assay = NULL, layer = "data") {
  version <- detect_seurat_version()
  
  if (is.null(assay)) {
    if (version >= 5) {
      assay <- SeuratObject::DefaultAssay(object)
    } else {
      assay <- Seurat::DefaultAssay(object)
    }
  }
  
  if (version >= 5) {
    # Seurat v5 uses layers
    tryCatch({
      data <- SeuratObject::LayerData(object, assay = assay, layer = layer)
    }, error = function(e) {
      # Fallback for v5 objects that might use old slots
      data <- SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    })
  } else {
    # Seurat v4 uses slots
    data <- Seurat::GetAssayData(object, assay = assay, slot = layer)
  }
  
  return(data)
}

#' Set Data in Seurat Object
#' @description Sets data in Seurat object, compatible with v4 and v5.
#' @param object Seurat object
#' @param data Matrix to set
#' @param assay Assay name
#' @param layer Layer/slot name
#' @return Modified Seurat object
#' @keywords internal
set_seurat_data <- function(object, data, assay, layer = "data") {
  version <- detect_seurat_version()
  
  if (version >= 5) {
    # Seurat v5
    tryCatch({
      SeuratObject::LayerData(object, assay = assay, layer = layer) <- data
    }, error = function(e) {
      SeuratObject::SetAssayData(object, assay = assay, slot = layer, new.data = data)
    })
  } else {
    # Seurat v4
    object <- Seurat::SetAssayData(object, assay = assay, slot = layer, new.data = data)
  }
  
  return(object)
}

#' MAGIC for Seurat Objects
#'
#' @description Run MAGIC on a Seurat object and store results in a new assay.
#'
#' @param data A Seurat object.
#' @param assay Character. Name of the assay to use as input. Default is the
#'   default assay of the object.
#' @param slot Character. Slot to use for input data. Default is "data"
#'   (normalized data). For Seurat v5, this corresponds to the layer name.
#' @param genes Character vector of gene names to return, or "all_genes" to
#'   return all genes. Default is "all_genes".
#' @param new_assay_name Character. Name for the new assay containing MAGIC
#'   results. Default is "MAGIC_{assay}".
#' @param knn Integer. Number of nearest neighbors. Default is 5.
#' @param knn_max Integer or NULL. Maximum neighbors. Default is NULL.
#' @param decay Numeric. Alpha decay parameter. Default is 1.
#' @param t Integer or "auto". Diffusion steps. Default is 3.
#' @param npca Integer. Number of PCA components. Default is 100.
#' @param solver Character. "exact" or "approximate". Default is "exact".
#' @param init Previous MAGIC result. Default is NULL.
#' @param t_max Integer. Maximum t for auto selection. Default is 20.
#' @param distance Character. Distance metric. Default is "euclidean".
#' @param n_jobs Integer. Number of parallel workers. Default is 1.
#' @param seed Integer or NULL. Random seed. Default is NULL.
#' @param verbose Logical or integer. Verbosity. Default is TRUE.
#' @param ... Additional arguments passed to magic.matrix.
#'
#' @return The Seurat object with a new assay containing MAGIC-imputed values.
#'   The new assay is named "MAGIC_{original_assay}" by default.
#'
#' @details This function:
#'   \enumerate{
#'     \item Extracts normalized data from the specified assay
#'     \item Transposes to cells x genes format
#'     \item Runs MAGIC imputation
#'     \item Creates a new assay with the imputed data
#'     \item Stores MAGIC parameters in the object's misc slot
#'   }
#'
#'   The function is compatible with both Seurat v4 and v5. For v4, it uses
#'   GetAssayData/SetAssayData. For v5, it uses LayerData where available.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Create example Seurat object
#' counts <- matrix(rpois(3000, lambda = 5), nrow = 30, ncol = 100)
#' rownames(counts) <- paste0("Gene", 1:30)
#' colnames(counts) <- paste0("Cell", 1:100)
#'
#' seurat_obj <- CreateSeuratObject(counts = counts)
#' seurat_obj <- NormalizeData(seurat_obj)
#'
#' # Run MAGIC
#' seurat_obj <- magic(seurat_obj, t = 3)
#'
#' # Access MAGIC results
#' GetAssayData(seurat_obj, assay = "MAGIC_RNA", slot = "data")
#'
#' # Make MAGIC assay the default
#' DefaultAssay(seurat_obj) <- "MAGIC_RNA"
#' }
#'
#' @rdname magic
#' @export
magic.Seurat <- function(data,
                          assay = NULL,
                          slot = "data",
                          genes = "all_genes",
                          new_assay_name = NULL,
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
  
  solver <- match.arg(solver)
  distance <- match.arg(distance)
  
  # Check Seurat is available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for this method.", call. = FALSE)
  }
  
  version <- detect_seurat_version()
  vmessage(sprintf("Detected Seurat v%d", version), verbose)
  
  # Get default assay if not specified
  if (is.null(assay)) {
    if (version >= 5) {
      assay <- SeuratObject::DefaultAssay(data)
    } else {
      assay <- Seurat::DefaultAssay(data)
    }
  }
  
  # Set default new assay name
  if (is.null(new_assay_name)) {
    new_assay_name <- paste0("MAGIC_", assay)
  }
  
  vmessage(sprintf("Running MAGIC on Seurat object (assay: %s)...", assay), verbose)
  
  # Extract data (genes x cells)
  expr_data <- get_seurat_data(data, assay = assay, layer = slot)
  
  # Transpose to cells x genes for MAGIC
  expr_data_t <- Matrix::t(expr_data)
  
  # Convert to dense if sparse
  if (inherits(expr_data_t, "sparseMatrix")) {
    expr_data_t <- as.matrix(expr_data_t)
  }
  
  # Run MAGIC
  magic_result <- magic.matrix(
    data = expr_data_t,
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
  
  # Get imputed data and transpose back to genes x cells
  imputed_data <- t(as.matrix(magic_result$result))
  
  # Create new assay
  vmessage(sprintf("Creating new assay '%s'...", new_assay_name), verbose)
  
  if (version >= 5) {
    # Seurat v5
    tryCatch({
      new_assay <- SeuratObject::CreateAssay5Object(data = imputed_data)
    }, error = function(e) {
      # Fallback to v4 style
      new_assay <- SeuratObject::CreateAssayObject(data = imputed_data)
    })
  } else {
    # Seurat v4
    new_assay <- Seurat::CreateAssayObject(data = imputed_data)
  }
  
  # Add assay to object
  data[[new_assay_name]] <- new_assay
  
  # Store MAGIC info in misc slot
  if (is.null(data@misc$MAGIC)) {
    data@misc$MAGIC <- list()
  }
  data@misc$MAGIC[[new_assay_name]] <- list(
    params = magic_result$params,
    run_time = Sys.time()
  )
  
  vmessage(sprintf("MAGIC complete. Results stored in assay '%s'.", new_assay_name),
           verbose)
  vmessage(sprintf("To use MAGIC results, run: DefaultAssay(object) <- '%s'",
                   new_assay_name), verbose)
  
  return(data)
}

#' MAGIC for SingleCellExperiment Objects
#'
#' @description Run MAGIC on a SingleCellExperiment object.
#'
#' @param data A SingleCellExperiment object.
#' @param assay_name Character. Name of the assay to use. Default is "logcounts".
#' @param genes Gene selection. Default is "all_genes".
#' @param new_assay_name Name for the new assay. Default is "MAGIC".
#' @param knn Number of neighbors. Default is 5.
#' @param knn_max Maximum neighbors. Default is NULL.
#' @param decay Alpha decay. Default is 1.
#' @param t Diffusion steps. Default is 3.
#' @param npca PCA components. Default is 100.
#' @param solver Solver type. Default is "exact".
#' @param init Previous result. Default is NULL.
#' @param t_max Maximum t. Default is 20.
#' @param distance Distance metric. Default is "euclidean".
#' @param n_jobs Parallel workers. Default is 1.
#' @param seed Random seed. Default is NULL.
#' @param verbose Verbosity. Default is TRUE.
#' @param ... Additional arguments.
#'
#' @return SingleCellExperiment object with MAGIC results in a new assay.
#'
#' @rdname magic
#' @export
magic.SingleCellExperiment <- function(data,
                                        assay_name = "logcounts",
                                        genes = "all_genes",
                                        new_assay_name = "MAGIC",
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
  
  solver <- match.arg(solver)
  distance <- match.arg(distance)
  
  # Check packages
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package is required.", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required.", call. = FALSE)
  }
  
  vmessage("Running MAGIC on SingleCellExperiment object...", verbose)
  
  # Get assay data (genes x cells)
  if (!assay_name %in% SummarizedExperiment::assayNames(data)) {
    stop(sprintf("Assay '%s' not found. Available: %s",
                 assay_name,
                 paste(SummarizedExperiment::assayNames(data), collapse = ", ")),
         call. = FALSE)
  }
  
  expr_data <- SummarizedExperiment::assay(data, assay_name)
  
  # Transpose to cells x genes
  expr_data_t <- Matrix::t(expr_data)
  if (inherits(expr_data_t, "sparseMatrix")) {
    expr_data_t <- as.matrix(expr_data_t)
  }
  
  # Run MAGIC
  magic_result <- magic.matrix(
    data = expr_data_t,
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
  
  # Transpose back and store
  imputed_data <- t(as.matrix(magic_result$result))
  
  # Add as new assay
  SummarizedExperiment::assay(data, new_assay_name) <- imputed_data
  
  # Store parameters in metadata
  S4Vectors::metadata(data)$MAGIC <- list(
    params = magic_result$params,
    run_time = Sys.time()
  )
  
  vmessage(sprintf("MAGIC complete. Results stored in assay '%s'.", new_assay_name),
           verbose)
  

  return(data)
}

#' Extract MAGIC Results from Seurat Object
#'
#' @description Convenience function to extract MAGIC-imputed data from a
#'   Seurat object.
#'
#' @param object Seurat object that has been processed with magic()
#' @param assay Name of the MAGIC assay. Default is "MAGIC_RNA".
#' @param genes Genes to extract. Default is NULL (all genes).
#' @param cells Cells to extract. Default is NULL (all cells).
#' @param as_matrix Logical. Return as matrix instead of data.frame? Default TRUE.
#'
#' @return Matrix or data.frame of MAGIC-imputed expression values (cells x genes).
#'
#' @export
#' @examples
#' \dontrun{
#' # After running magic() on Seurat object
#' imputed_data <- GetMAGICData(seurat_obj)
#' imputed_subset <- GetMAGICData(seurat_obj, genes = c("Gene1", "Gene2"))
#' }
GetMAGICData <- function(object,
                          assay = "MAGIC_RNA",
                          genes = NULL,
                          cells = NULL,
                          as_matrix = TRUE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.", call. = FALSE)
  }
  
  if (!assay %in% names(object@assays)) {
    stop(sprintf("Assay '%s' not found. Available MAGIC assays: %s",
                 assay,
                 paste(grep("^MAGIC", names(object@assays), value = TRUE),
                       collapse = ", ")),
         call. = FALSE)
  }
  
  # Get data
  data <- get_seurat_data(object, assay = assay, layer = "data")
  
  # Subset genes
  if (!is.null(genes)) {
    genes <- genes[genes %in% rownames(data)]
    if (length(genes) == 0) {
      stop("No matching genes found.", call. = FALSE)
    }
    data <- data[genes, , drop = FALSE]
  }
  
  # Subset cells
  if (!is.null(cells)) {
    cells <- cells[cells %in% colnames(data)]
    if (length(cells) == 0) {
      stop("No matching cells found.", call. = FALSE)
    }
    data <- data[, cells, drop = FALSE]
  }
  
  # Transpose to cells x genes and convert
  data <- Matrix::t(data)
  
  if (as_matrix) {
    return(as.matrix(data))
  } else {
    return(as.data.frame(as.matrix(data)))
  }
}
