# MAGICR

<!-- badges: start -->
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/MAGICR)](https://zaoqu-liu.r-universe.dev/MAGICR)
[![R-CMD-check](https://github.com/Zaoqu-Liu/MAGICR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zaoqu-Liu/MAGICR/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/MAGICR)](https://CRAN.R-project.org/package=MAGICR)
[![License: GPL-2](https://img.shields.io/badge/License-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
<!-- badges: end -->

📖 **Documentation**: <https://zaoqu-liu.github.io/MAGICR>

## Overview

**MAGICR** is a native R implementation of the MAGIC (**M**arkov **A**ffinity-based **G**raph **I**mputation of **C**ells) algorithm for denoising and imputation of single-cell RNA sequencing (scRNA-seq) data.

Single-cell RNA sequencing enables transcriptomic profiling at single-cell resolution, but technical limitations result in sparse count matrices with prevalent dropout events. MAGIC addresses this challenge by leveraging the underlying manifold structure of the data through diffusion geometry, effectively recovering missing transcripts while preserving biological signal.

## Algorithm

MAGIC performs data imputation through the following steps:

1. **Dimensionality Reduction**: Optional PCA to reduce computational complexity
2. **Graph Construction**: Build a k-nearest neighbor (kNN) graph based on cell-cell distances
3. **Kernel Computation**: Construct an α-decaying kernel with adaptive bandwidth:

$$K(x_i, x_j) = \exp\left(-\left(\frac{\|x_i - x_j\|_2}{\varepsilon_i}\right)^\alpha\right)$$

where $\varepsilon_i$ is the distance to the $k$-th nearest neighbor of cell $i$.

4. **Diffusion Operator**: Compute the row-stochastic Markov transition matrix:

$$P = D^{-1}K$$

5. **Imputation**: Apply powered diffusion:

$$X_{\text{imputed}} = P^t X$$

The diffusion time $t$ controls the degree of smoothing. When `t = "auto"`, optimal $t$ is selected by minimizing Procrustes disparity between successive iterations.

## Installation

### From R-universe (Recommended)

```r
install.packages("MAGICR", repos = "https://zaoqu-liu.r-universe.dev")
```

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("Zaoqu-Liu/MAGICR")
```

### Dependencies

**Required:**
- R (≥ 3.6.0)
- Matrix, irlba, RANN, future

**Optional:**
- Seurat (≥ 4.0.0) for Seurat object integration
- SingleCellExperiment for Bioconductor compatibility
- RcppArmadillo for C++ acceleration

## Usage

### Basic Usage

```r
library(MAGICR)

# Load example data
data(magic_testdata)

# Run MAGIC with default parameters
result <- magic(magic_testdata, t = 3)

# Access imputed data
imputed_data <- as.matrix(result)
```

### With Seurat Objects

```r
library(Seurat)
library(MAGICR)

# Run MAGIC on normalized data
seurat_obj <- magic(seurat_obj, assay = "RNA", t = 3)

# Results stored in "MAGIC_RNA" assay
DefaultAssay(seurat_obj) <- "MAGIC_RNA"
```

### Automatic Diffusion Time Selection

```r
# Automatically determine optimal t using Procrustes analysis
result <- magic(data, t = "auto", t_max = 20)

# Check selected t
print(result$params$t)
```

### Parameter Tuning

| Parameter | Default | Description |
|-----------|---------|-------------|
| `knn` | 5 | Number of nearest neighbors for bandwidth estimation |
| `knn_max` | 15 | Maximum neighbors for graph construction |
| `decay` | 1 | α parameter controlling kernel sharpness |
| `t` | 3 | Diffusion time (integer or "auto") |
| `npca` | 100 | PCA components for distance computation |
| `solver` | "exact" | "exact" or "approximate" (PCA-space) |
| `distance` | "euclidean" | Distance metric: "euclidean", "cosine", "correlation" |

## Methodological Details

### Adaptive Bandwidth

The α-decaying kernel uses cell-specific bandwidth $\varepsilon_i$, defined as the distance to the $k$-th nearest neighbor. This provides:
- Local density adaptation
- Robustness to varying cell densities
- Preservation of rare cell populations

### Diffusion Time Selection

When `t = "auto"`, MAGIC iteratively applies the diffusion operator and monitors convergence using Procrustes disparity:

$$d_{\text{Procrustes}}(X_t, X_{t+1}) < \theta$$

The algorithm terminates when the disparity falls below threshold $\theta$ (default: 0.001).

### Computational Considerations

- **Exact solver**: Computes imputation in full gene space. Recommended for smaller datasets or when preserving all gene relationships is critical.
- **Approximate solver**: Operates in PCA space and projects back. Faster for large datasets but may introduce minor approximation errors.

## Citation

If you use MAGICR in your research, please cite:

> van Dijk, D., Sharma, R., Nainys, J., Yim, K., Kathail, P., Carr, A.J., Burdziak, C., Moon, K.R., Chaffer, C.L., Pattabiraman, D., Bierie, B., Mazutis, L., Wolf, G., Krishnaswamy, S., & Pe'er, D. (2018). **Recovering Gene Interactions from Single-Cell Data Using Data Diffusion**. *Cell*, 174(3), 716-729.e27. https://doi.org/10.1016/j.cell.2018.05.061

## Related Resources

- [Original Python MAGIC](https://github.com/KrishnaswamyLab/MAGIC) - Krishnaswamy Lab implementation
- [PHATE](https://github.com/KrishnaswamyLab/PHATE) - Dimensionality reduction using similar diffusion principles
- [graphtools](https://github.com/KrishnaswamyLab/graphtools) - Python library for graph-based data analysis

## License

GPL-2

## Author

**Zaoqu Liu** (liuzaoqu@163.com)

- GitHub: [@Zaoqu-Liu](https://github.com/Zaoqu-Liu)
- ORCID: [0000-0002-0452-742X](https://orcid.org/0000-0002-0452-742X)

## Acknowledgments

This implementation is based on the original MAGIC algorithm developed by the [Krishnaswamy Lab](https://www.krishnaswamylab.org/) at Yale University. The MATLAB implementation served as the primary reference for this native R port.
