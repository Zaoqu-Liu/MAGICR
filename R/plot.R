#' @title MAGIC Plotting Functions
#' @description Visualization functions for MAGIC results.
#' @name plot
#' @importFrom grDevices adjustcolor colorRamp hcl.colors rgb
#' @importFrom graphics hist par
#' @importFrom stats var
NULL

#' Animate MAGIC Diffusion
#'
#' @description Create an animation showing how gene-gene relationships
#'   evolve with increasing diffusion time t.
#'
#' @param data Input data matrix (cells x genes).
#' @param gene_x Gene for x-axis (name or index).
#' @param gene_y Gene for y-axis (name or index).
#' @param gene_color Gene for coloring points (optional).
#' @param t_max Maximum t value. Default is 20.
#' @param delay Number of initial frames before diffusion. Default is 2.
#' @param operator Pre-computed magic object. Default is NULL.
#' @param filename Output filename (.gif). Default is NULL (no save).
#' @param width Plot width in pixels. Default is 600.
#' @param height Plot height in pixels. Default is 500.
#' @param point_size Point size. Default is 1.
#' @param interval Milliseconds between frames. Default is 200.
#' @param verbose Verbosity. Default is FALSE.
#' @param ... Additional arguments passed to magic().
#'
#' @return A list of ggplot objects (one per frame), or writes a GIF file.
#'
#' @details This function creates an animation showing how MAGIC diffusion
#'   progressively reveals gene-gene relationships. It requires the 
#'   \code{gganimate} and \code{gifski} packages for GIF output.
#'
#' @export
#' @examples
#' \dontrun{
#' data(magic_testdata)
#' 
#' # Create animation frames
#' frames <- animate_magic(magic_testdata, 
#'                         gene_x = 1, gene_y = 2,
#'                         t_max = 10)
#' 
#' # Save as GIF (requires gganimate and gifski)
#' animate_magic(magic_testdata, 
#'               gene_x = 1, gene_y = 2,
#'               t_max = 10,
#'               filename = "magic_animation.gif")
#' }
animate_magic <- function(data,
                          gene_x,
                          gene_y,
                          gene_color = NULL,
                          t_max = 20,
                          delay = 2,
                          operator = NULL,
                          filename = NULL,
                          width = 600,
                          height = 500,
                          point_size = 1,
                          interval = 200,
                          verbose = FALSE,
                          ...) {
  
  # Validate inputs
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a matrix or data.frame", call. = FALSE)
  }
  
  data <- as.matrix(data)
  
  # Resolve gene names/indices
  gene_x <- resolve_gene(gene_x, colnames(data))
  gene_y <- resolve_gene(gene_y, colnames(data))
  if (!is.null(gene_color)) {
    gene_color <- resolve_gene(gene_color, colnames(data))
  }
  
  # Build or use operator
  if (is.null(operator)) {
    vmessage("Building MAGIC operator...", verbose)
    operator <- magic(data, t = 1, verbose = verbose, ...)
  }
  
  # Generate frames
  frames <- list()
  
  # Original data (t = 0)
  for (i in seq_len(delay)) {
    frames[[length(frames) + 1]] <- list(
      t = 0,
      x = data[, gene_x],
      y = data[, gene_y],
      color = if (!is.null(gene_color)) data[, gene_color] else NULL
    )
  }
  
  # Diffused data
  for (t in seq_len(t_max)) {
    vmessage(sprintf("  Computing t = %d...", t), verbose)
    
    # Re-run magic with this t
    result <- magic(data, t = t, init = operator, verbose = FALSE, ...)
    imputed <- as.matrix(result)
    
    frames[[length(frames) + 1]] <- list(
      t = t,
      x = imputed[, gene_x],
      y = imputed[, gene_y],
      color = if (!is.null(gene_color)) imputed[, gene_color] else NULL
    )
  }
  
  # Get gene names for labels
  gene_x_name <- if (is.character(gene_x)) gene_x else colnames(data)[gene_x]
  gene_y_name <- if (is.character(gene_y)) gene_y else colnames(data)[gene_y]
  gene_color_name <- if (!is.null(gene_color)) {
    if (is.character(gene_color)) gene_color else colnames(data)[gene_color]
  } else NULL
  
  # Create plots
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- create_animation_plots(frames, gene_x_name, gene_y_name, 
                                    gene_color_name, point_size)
    
    # Save GIF if filename provided
    if (!is.null(filename)) {
      save_animation_gif(plots, filename, width, height, interval)
    }
    
    return(invisible(plots))
  } else {
    # Base R fallback - just create a function to plot each frame
    plot_frame <- function(frame_idx) {
      frame <- frames[[frame_idx]]
      if (!is.null(frame$color)) {
        plot(frame$x, frame$y, 
             col = color_scale(frame$color),
             pch = 16, cex = point_size,
             main = sprintf("t = %d", frame$t),
             xlab = gene_x_name, ylab = gene_y_name)
      } else {
        plot(frame$x, frame$y,
             col = adjustcolor("steelblue", 0.5),
             pch = 16, cex = point_size,
             main = sprintf("t = %d", frame$t),
             xlab = gene_x_name, ylab = gene_y_name)
      }
    }
    
    return(list(
      frames = frames,
      plot_frame = plot_frame,
      n_frames = length(frames)
    ))
  }
}

#' Resolve gene name or index
#' @keywords internal
resolve_gene <- function(gene, colnames) {
  if (is.character(gene)) {
    if (!gene %in% colnames) {
      stop(sprintf("Gene '%s' not found", gene), call. = FALSE)
    }
    return(gene)
  } else if (is.numeric(gene)) {
    gene <- as.integer(gene)
    if (gene < 1 || gene > length(colnames)) {
      stop(sprintf("Gene index %d out of range", gene), call. = FALSE)
    }
    return(gene)
  }
  stop("Gene must be a name or index", call. = FALSE)
}

#' Create animation plots with ggplot2
#' @keywords internal
create_animation_plots <- function(frames, gene_x_name, gene_y_name, 
                                   gene_color_name, point_size) {
  # Avoid R CMD check note
  x <- y <- color <- NULL
  
  plots <- list()
  
  for (i in seq_along(frames)) {
    frame <- frames[[i]]
    df <- data.frame(x = frame$x, y = frame$y)
    
    if (!is.null(frame$color)) {
      df$color <- frame$color
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = color)) +
        ggplot2::geom_point(size = point_size, alpha = 0.6) +
        ggplot2::scale_color_viridis_c(name = gene_color_name)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(size = point_size, alpha = 0.6, color = "steelblue")
    }
    
    p <- p +
      ggplot2::labs(
        title = sprintf("MAGIC: t = %d", frame$t),
        x = gene_x_name,
        y = gene_y_name
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.title = ggplot2::element_text(size = 12)
      )
    
    plots[[i]] <- p
  }
  
  return(plots)
}

#' Save animation as GIF
#' @keywords internal
save_animation_gif <- function(plots, filename, width, height, interval) {
  if (!requireNamespace("gifski", quietly = TRUE)) {
    warning("Package 'gifski' required to save GIF. Install with: ",
            "install.packages('gifski')")
    return(invisible(NULL))
  }
  
  # Create temporary PNG files
  temp_dir <- tempdir()
  png_files <- character(length(plots))
  
  for (i in seq_along(plots)) {
    png_files[i] <- file.path(temp_dir, sprintf("frame_%04d.png", i))
    ggplot2::ggsave(png_files[i], plots[[i]], 
                    width = width / 100, height = height / 100, dpi = 100)
  }
  
  # Create GIF
  gifski::gifski(png_files, gif_file = filename, 
                 width = width, height = height,
                 delay = interval / 1000)
  
  # Clean up
 unlink(png_files)
  
  message(sprintf("Animation saved to: %s", filename))
  invisible(filename)
}

#' Simple color scale
#' @keywords internal
color_scale <- function(values) {
  values <- (values - min(values)) / (max(values) - min(values) + .Machine$double.eps)
  rgb(colorRamp(c("yellow", "red"))(values) / 255)
}

#' Plot Gene-Gene Scatter
#'
#' @description Create a scatter plot of two genes before and after MAGIC.
#'
#' @param data Original data matrix.
#' @param magic_result MAGIC result object.
#' @param gene_x Gene for x-axis.
#' @param gene_y Gene for y-axis.
#' @param gene_color Gene for coloring (optional).
#' @param point_size Point size. Default is 1.
#'
#' @return A ggplot object or plots using base R.
#'
#' @export
#' @examples
#' \dontrun{
#' data(magic_testdata)
#' result <- magic(magic_testdata, t = 3)
#' plot_magic_genes(magic_testdata, result, 1, 2)
#' }
plot_magic_genes <- function(data, magic_result, gene_x, gene_y, 
                              gene_color = NULL, point_size = 1) {
  # Avoid R CMD check note
  x <- y <- type <- color <- density <- NULL
  
  data <- as.matrix(data)
  imputed <- as.matrix(magic_result)
  
  gene_x <- resolve_gene(gene_x, colnames(data))
  gene_y <- resolve_gene(gene_y, colnames(data))
  
  gene_x_name <- if (is.character(gene_x)) gene_x else colnames(data)[gene_x]
  gene_y_name <- if (is.character(gene_y)) gene_y else colnames(data)[gene_y]
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df_orig <- data.frame(
      x = data[, gene_x],
      y = data[, gene_y],
      type = "Original"
    )
    df_magic <- data.frame(
      x = imputed[, gene_x],
      y = imputed[, gene_y],
      type = "After MAGIC"
    )
    df <- rbind(df_orig, df_magic)
    df$type <- factor(df$type, levels = c("Original", "After MAGIC"))
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(size = point_size, alpha = 0.5, color = "steelblue") +
      ggplot2::facet_wrap(~ type, scales = "free") +
      ggplot2::labs(x = gene_x_name, y = gene_y_name) +
      ggplot2::theme_minimal()
    
    print(p)
    return(invisible(p))
  } else {
    par(mfrow = c(1, 2))
    
    plot(data[, gene_x], data[, gene_y],
         pch = 16, col = adjustcolor("steelblue", 0.5),
         cex = point_size, main = "Original",
         xlab = gene_x_name, ylab = gene_y_name)
    
    plot(imputed[, gene_x], imputed[, gene_y],
         pch = 16, col = adjustcolor("steelblue", 0.5),
         cex = point_size, main = "After MAGIC",
         xlab = gene_x_name, ylab = gene_y_name)
    
    par(mfrow = c(1, 1))
  }
}
