#' Example Single-Cell Data for MAGICR
#'
#' @description A small example dataset for demonstrating MAGICR functionality.
#'   Contains 500 cells and 100 genes with simulated count data.
#'
#' @format A matrix with 500 rows (cells) and 100 columns (genes):
#' \describe{
#'   \item{rows}{Cell identifiers (Cell_1 to Cell_500)}
#'   \item{columns}{Gene identifiers (Gene_1 to Gene_100)}
#'   \item{values}{Integer count values (Poisson distributed with lambda=5)}
#' }
#'
#' @source Simulated data for package demonstration
#'
#' @examples
#' \dontrun{
#' data(magic_testdata)
#' dim(magic_testdata)
#' # [1] 500 100
#'
#' # Run MAGIC on example data
#' result <- magic(magic_testdata, t = 3)
#' }
"magic_testdata"
