## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.align = "center",
  message = FALSE,
  warning = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
# # From R-universe (recommended)
# install.packages("MAGICR", repos = "https://zaoqu-liu.r-universe.dev")
# 
# # From GitHub
# remotes::install_github("Zaoqu-Liu/MAGICR")

## -----------------------------------------------------------------------------
library(MAGICR)

# Load example data
data(magic_testdata)

# Check data dimensions
dim(magic_testdata)

# Preview the data
magic_testdata[1:5, 1:5]

## -----------------------------------------------------------------------------
# Run MAGIC with default parameters
result <- magic(magic_testdata, t = 3)

# View result summary
print(result)

## -----------------------------------------------------------------------------
# Get imputed matrix
imputed_data <- as.matrix(result)

# Compare original vs imputed
cat("Original data range:", range(magic_testdata), "\n")
cat("Imputed data range:", range(imputed_data), "\n")

## ----fig.width=8, fig.height=4------------------------------------------------
par(mfrow = c(1, 2))

# Original data distribution
hist(as.vector(as.matrix(magic_testdata)), 
     breaks = 50, 
     main = "Original Data Distribution",
     xlab = "Expression", 
     col = "#3498db", 
     border = "white")

# Imputed data distribution
hist(as.vector(imputed_data), 
     breaks = 50, 
     main = "Imputed Data Distribution",
     xlab = "Expression", 
     col = "#e74c3c", 
     border = "white")

## ----fig.width=6, fig.height=6------------------------------------------------
# Select two genes
gene1 <- colnames(magic_testdata)[1]
gene2 <- colnames(magic_testdata)[2]

par(mfrow = c(1, 2))

# Original
plot(magic_testdata[, gene1], magic_testdata[, gene2],
     pch = 16, col = adjustcolor("#3498db", 0.5),
     xlab = gene1, ylab = gene2,
     main = "Original")

# Imputed
plot(imputed_data[, gene1], imputed_data[, gene2],
     pch = 16, col = adjustcolor("#e74c3c", 0.5),
     xlab = gene1, ylab = gene2,
     main = "After MAGIC")

## -----------------------------------------------------------------------------
# Automatic t selection
result_auto <- magic(magic_testdata, t = "auto", t_max = 10)
cat("Automatically selected t:", result_auto$params$t, "\n")

## -----------------------------------------------------------------------------
sessionInfo()

