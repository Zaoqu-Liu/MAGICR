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

## ----fig.width=8, fig.height=4------------------------------------------------
# Visualize α-decaying kernel
x <- seq(0, 3, length.out = 100)

par(mfrow = c(1, 2))

# Effect of distance
plot(x, exp(-x), type = "l", lwd = 2, col = "#e74c3c",
     xlab = "Normalized Distance (d/ε)", ylab = "Affinity K(x,y)",
     main = "α-Decaying Kernel (α = 1)")
abline(h = 0.5, lty = 2, col = "gray")
abline(v = log(2), lty = 2, col = "gray")

# Effect of α
plot(x, exp(-x^0.5), type = "l", lwd = 2, col = "#3498db",
     xlab = "Normalized Distance (d/ε)", ylab = "Affinity K(x,y)",
     main = "Effect of α Parameter", ylim = c(0, 1))
lines(x, exp(-x^1), lwd = 2, col = "#e74c3c")
lines(x, exp(-x^2), lwd = 2, col = "#2ecc71")
legend("topright", legend = c("α = 0.5", "α = 1", "α = 2"),
       col = c("#3498db", "#e74c3c", "#2ecc71"), lwd = 2)

## ----fig.width=8, fig.height=4------------------------------------------------
# Simulate diffusion effect
set.seed(42)
n <- 100
x <- c(rnorm(50, -2, 0.5), rnorm(50, 2, 0.5))
y <- c(rnorm(50, 0, 0.5), rnorm(50, 0, 0.5))

# Add noise
x_noisy <- x + rnorm(n, 0, 0.3)
y_noisy <- y + rnorm(n, 0, 0.3)

par(mfrow = c(1, 2))
plot(x_noisy, y_noisy, pch = 16, col = adjustcolor("#3498db", 0.6),
     main = "Before Diffusion (Noisy)", xlab = "Gene 1", ylab = "Gene 2")

# Simple smoothing simulation
x_smooth <- 0.7 * x + 0.3 * x_noisy
y_smooth <- 0.7 * y + 0.3 * y_noisy
plot(x_smooth, y_smooth, pch = 16, col = adjustcolor("#e74c3c", 0.6),
     main = "After Diffusion (Denoised)", xlab = "Gene 1", ylab = "Gene 2")

## -----------------------------------------------------------------------------
library(MAGICR)
data(magic_testdata)

# Run with automatic t selection
result <- magic(magic_testdata, t = "auto", t_max = 15, verbose = TRUE)

cat("\nOptimal t selected:", result$params$t, "\n")

## -----------------------------------------------------------------------------
sessionInfo()

