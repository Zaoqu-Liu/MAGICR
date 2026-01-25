# MAGICR <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Zaoqu-Liu/MAGICR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zaoqu-Liu/MAGICR/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/MAGICR)](https://CRAN.R-project.org/package=MAGICR)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

**MAGICR** is a native R implementation of **MAGIC** (Markov Affinity-based Graph Imputation of Cells) for denoising and imputation of single-cell RNA sequencing data.

MAGIC uses diffusion geometry to denoise single-cell data and fill in missing transcripts, as described in [van Dijk et al., 2018](https://doi.org/10.1016/j.cell.2018.05.061).

## ✨ Features

- **Pure R implementation** - No Python dependencies required
- **C++ acceleration** - Optional RcppArmadillo backend for speed
- **Parallel computing** - Built-in support via the `future` framework
- **Seurat integration** - Seamless support for Seurat v4 and v5 objects
- **SingleCellExperiment support** - Works with Bioconductor workflows
- **Cross-platform** - Works on Windows, macOS, and Linux

## 📦 Installation

### From GitHub (recommended)

```r
# Install devtools if needed
install.packages("devtools")

# Install MAGICR
devtools::install_github("Zaoqu-Liu/MAGICR")
```

### Dependencies

Required:
- R (>= 3.6.0)
- Matrix
- irlba
- RANN
- future
- future.apply
- progressr

Optional (for Seurat integration):
- Seurat (>= 4.0.0)
- SeuratObject

Optional (for Bioconductor integration):
- SingleCellExperiment
- SummarizedExperiment

## 🚀 Quick Start

### Basic Usage

```r
library(MAGICR)

# Create example data (500 cells x 100 genes)
set.seed(42)
data <- matrix(rpois(50000, lambda = 2), nrow = 500, ncol = 100)
colnames(data) <- paste0("Gene", 1:100)
rownames(data) <- paste0("Cell", 1:500)

# Run MAGIC with default parameters
result <- magic(data, t = 3, verbose = TRUE)

# Access imputed data
imputed <- as.matrix(result)

# Or as data.frame
imputed_df <- as.data.frame(result)
```

### With Seurat Objects

```r
library(Seurat)
library(MAGICR)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = t(data))
seurat_obj <- NormalizeData(seurat_obj)

# Run MAGIC
seurat_obj <- magic(seurat_obj, t = 3)

# Results are stored in "MAGIC_RNA" assay
DefaultAssay(seurat_obj) <- "MAGIC_RNA"

# Visualize
VlnPlot(seurat_obj, features = c("Gene1", "Gene2"))
```

### With SingleCellExperiment

```r
library(SingleCellExperiment)
library(MAGICR)

# Run MAGIC on SCE object
sce <- magic(sce, assay_name = "logcounts", t = 3)

# Results in "MAGIC" assay
assay(sce, "MAGIC")
```

## 📖 Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `knn` | 5 | Number of nearest neighbors for kernel bandwidth |
| `knn_max` | NULL | Maximum neighbors (default: 3 × knn) |
| `decay` | 1 | Alpha decay parameter for kernel |
| `t` | 3 | Number of diffusion steps (or "auto") |
| `npca` | 100 | Number of PCA components for distance calculation |
| `solver` | "exact" | "exact" or "approximate" |
| `distance` | "euclidean" | Distance metric: "euclidean", "cosine", "manhattan", "correlation" |
| `n_jobs` | 1 | Number of parallel workers (-1 for all) |

## ⚡ Performance Tips

### Enable Parallel Computing

```r
library(future)

# Use all available cores
plan(multisession, workers = parallel::detectCores())

# Run MAGIC
result <- magic(data, t = 3, n_jobs = -1)
```

### For Large Datasets

```r
# Use approximate solver (faster, uses PCA space)
result <- magic(data, 
                t = 3, 
                solver = "approximate",
                npca = 50)
```

### Check C++ Acceleration

```r
# Check if C++ is available
has_cpp_acceleration()

# View configuration
magicr_config()
```

## 🔬 Algorithm Details

MAGIC performs the following steps:

1. **PCA Reduction** (optional): Reduces dimensions for faster distance computation
2. **kNN Graph Construction**: Builds a k-nearest neighbor graph
3. **Alpha-Decaying Kernel**: Computes adaptive bandwidth kernel
4. **Diffusion Operator**: Creates row-stochastic Markov matrix
5. **Diffusion**: Applies the operator t times to smooth expression

The key formula:

$$X_{imputed} = P^t \times X$$

where $P$ is the diffusion operator and $t$ controls the smoothing level.

### Automatic t Selection

When `t = "auto"`, MAGIC iteratively applies diffusion and stops when the Procrustes disparity between consecutive iterations falls below a threshold (default: 0.001).

## 📚 Citation

If you use MAGICR, please cite:

> van Dijk, D., Sharma, R., Nainys, J., et al. (2018). Recovering gene interactions from single-cell data using data diffusion. *Cell*, 174(3), 716-729.e27. https://doi.org/10.1016/j.cell.2018.05.061

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

GPL (>= 2)

## 👨‍💻 Author

**Zaoqu Liu** (liuzaoqu@163.com)

- GitHub: [@Zaoqu-Liu](https://github.com/Zaoqu-Liu)

## 🙏 Acknowledgments

- Original MAGIC algorithm by [Krishnaswamy Lab](https://github.com/KrishnaswamyLab/MAGIC)
- MATLAB implementation provided reference for this native R port
