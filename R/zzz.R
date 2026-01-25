#' @title Package Loading and Configuration
#' @description Internal package initialization functions.
#' @name zzz
#' @keywords internal
#' @import Matrix
#' @import methods
#' @importFrom stats median setNames
#' @importFrom graphics abline points plot
#' @importFrom utils packageVersion
#' @importFrom irlba prcomp_irlba
#' @importFrom RANN nn2
#' @importFrom Rcpp sourceCpp
#' @importFrom future plan nbrOfWorkers availableCores
#' @useDynLib MAGICR, .registration = TRUE
NULL

# Package environment for storing state
.magicr_env <- new.env(parent = emptyenv())

#' @importFrom utils packageVersion
.onLoad <- function(libname, pkgname) {
  # Check for optional C++ acceleration
  .magicr_env$use_rcpp <- tryCatch({
    .Call("_MAGICR_rcpp_available", PACKAGE = "MAGICR")
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  # Store package version
  .magicr_env$version <- utils::packageVersion("MAGICR")
  
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion("MAGICR")
  
  msg <- sprintf(
    "MAGICR version %s\n",
    version
  )
  
  # Check C++ availability
  if (isTRUE(.magicr_env$use_rcpp)) {
    msg <- paste0(msg, "  C++ acceleration: enabled\n")
  } else {
    msg <- paste0(msg, "  C++ acceleration: disabled (using pure R)\n")
  }
  
  # Check parallel capabilities
  if (requireNamespace("future", quietly = TRUE)) {
    msg <- paste0(msg, "  Parallel support: available (future package)\n")
  }
  
  packageStartupMessage(msg)
}

#' Check if C++ Acceleration is Available
#' 
#' @description Returns TRUE if C++ functions compiled successfully.
#' @return Logical
#' @export
#' @examples
#' has_cpp_acceleration()
has_cpp_acceleration <- function() {
  isTRUE(.magicr_env$use_rcpp)
}

#' Get MAGICR Configuration
#' 
#' @description Returns current configuration settings.
#' @return Named list of settings
#' @export
#' @examples
#' magicr_config()
magicr_config <- function() {
  list(
    version = as.character(.magicr_env$version %||% utils::packageVersion("MAGICR")),
    cpp_enabled = has_cpp_acceleration(),
    parallel_available = requireNamespace("future", quietly = TRUE),
    seurat_available = requireNamespace("Seurat", quietly = TRUE),
    sce_available = requireNamespace("SingleCellExperiment", quietly = TRUE)
  )
}

#' Set MAGICR Options
#' 
#' @description Configure package-wide options.
#' @param use_cpp Logical. Whether to use C++ acceleration when available.
#' @param verbose_default Integer. Default verbosity level (0, 1, or 2).
#' @return Invisibly returns previous settings
#' @export
#' @examples
#' \dontrun{
#' # Disable C++ acceleration
#' set_magicr_options(use_cpp = FALSE)
#' }
set_magicr_options <- function(use_cpp = NULL, verbose_default = NULL) {
  prev <- list(
    use_cpp = .magicr_env$use_rcpp,
    verbose_default = .magicr_env$verbose_default
  )
  
  if (!is.null(use_cpp)) {
    .magicr_env$use_rcpp <- use_cpp
  }
  
  if (!is.null(verbose_default)) {
    .magicr_env$verbose_default <- verbose_default
  }
  
  invisible(prev)
}
