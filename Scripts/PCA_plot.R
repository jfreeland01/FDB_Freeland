#' PCA plot
#'
#' Plots PCA from scores file (output of PCA_from_file)
#'
#' @param file File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param title Title of the plot
#' @param labels Show point labels (default = TRUE)
#' @param PCx, PCy Principal components to display
#' @param ellipse Draw confidence ellipses (default = FALSE)
#' @param conf Confidence level for ellipses (default = 0.95)
#' @param density Add marginal density plots (default = FALSE)
#' @param fliph, flipv Flip axes
#' @param show.legend Show legends (default = TRUE)
#' @param colors Named vector of colors per group
#' @param info.shape Vector of shape categories (same length/order as info.name)
#' @param shapes Named vector of shapes (e.g., c("V1" = 15, "V2" = 16))
#'
#' @export

library(vegan)
library(ggplot2)

PCA_plot <- function(file,
                     info.name,
                     info.type,
                     title = "",
                     labels = TRUE,
                     PCx = "PC1",
                     PCy = "PC2",
                     ellipse = FALSE,
                     conf = 0.95,
                     density = FALSE,
                     fliph = FALSE,
                     flipv = FALSE,
                     show.legend = TRUE,
                     colors = NULL,
                     info.shape = NULL,
                     shapes = NULL) {
  
  table <- read.table(file, header = TRUE)
  table$type <- info.type[match(table$Score, info.name)]
  
  if (!is.null(info.shape)) {
    table$shape <- info.shape[match(table$Score, info.name)]
  }
  
  if (grepl("scores_VARIMAX.txt", file)) {
    PCx <- gsub("PC", "V", PCx)
    PCy <- gsub("PC", "V", PCy)
  }
  
  if (fliph) table[[PCx]] <- -table[[PCx]]
  if (flipv) table[[PCy]] <- -table[[PCy]]
  
  sdev_name <- paste0(gsub("scores.txt", "", file), "sdev.txt")
  if (file.exists(sdev_name)) {
    sdev <- read.delim(sdev_name)$x
    pve <- round((sdev^2) / sum(sdev^2) * 100, 2)
    names(pve) <- paste0("PC", seq_along(pve))
    xlab <- if (!is.na(pve[PCx])) paste0(PCx, " (", pve[PCx], "%)") else PCx
    ylab <- if (!is.na(pve[PCy])) paste0(PCy, " (", pve[PCy], "%)") else PCy
  } else {
    xlab <- PCx
    ylab <- PCy
  }
  
  make_base_plot <- function(data, xvar, yvar, color_by, shape_by, title, xlab, ylab, show_legend, add_labels, colors, shapes) {
    if (!is.null(shape_by)) {
      p <- ggplot(data, aes_string(x = xvar, y = yvar, text = "Score")) +
        geom_point(aes(color = factor(color_by), shape = factor(shape_by)), size = 3, show.legend = show_legend)
      if (!is.null(shapes)) {
        p <- p + scale_shape_manual(values = shapes)
      }
      p <- p + guides(shape = guide_legend(title = "Batch"))
    } else {
      p <- ggplot(data, aes_string(x = xvar, y = yvar, text = "Score")) +
        geom_point(aes(color = factor(color_by)), shape = 16, size = 3, show.legend = show_legend)
    }
    
    if (!is.null(colors)) {
      p <- p + scale_color_manual(values = colors)
    }
    
    p <- p +
      labs(title = title, x = xlab, y = ylab) +
      guides(color = guide_legend(title = "Group")) +
      theme_bw(base_size = 18) +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.background = element_rect(),
        axis.text.x = element_text(margin = margin(b = 5)),
        axis.text.y = element_text(margin = margin(l = 5))
      )
    
    if (add_labels) {
      p <- p + geom_text(aes(label = Score), check_overlap = TRUE, size = 3)
    }
    
    return(p)
  }
  
  pcx.y <- make_base_plot(
    table, PCx, PCy,
    color_by = table$type,
    shape_by = if (!is.null(info.shape)) table$shape else NULL,
    title = title,
    xlab = xlab,
    ylab = ylab,
    show_legend = show.legend,
    add_labels = labels,
    colors = colors,
    shapes = shapes
  )
  
  if (ellipse) {
    plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, main = "", col = "white")
    ord <- vegan::ordiellipse(table[, c(PCx, PCy)], table$type, kind = "sd", conf = conf)
    
    cov_ellipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
      theta <- (0:npoints) * 2 * pi / npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }
    
    df_ell <- data.frame()
    for (g in levels(droplevels(factor(table$type)))) {
      coords <- cov_ellipse(ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale)
      colnames(coords) <- c(PCx, PCy)
      df_ell <- rbind(df_ell, cbind(as.data.frame(coords), type = g))
    }
    
    pcx.y <- pcx.y + geom_path(
      data = df_ell,
      aes_string(x = PCx, y = PCy, colour = "type"),
      size = 1,
      linetype = 1,
      inherit.aes = FALSE
    )
  }
  
  return(pcx.y)
}
