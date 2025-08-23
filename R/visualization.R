#' Plot mean and range of clusters for each particle
#'
#' Creates a scatter plot showing the mean and range (min/max) of data within clusters for each particle, colored by cluster.
#'
#' @param output_wvi A list containing at least elements \code{particles} (matrix) from WASABI output.
#' @param Y Numeric vector of data.
#' @return A ggplot2 plot object.
#' @export
#' @examples
#' \dontrun{
#' ggrange(output_wvi, Y)
#' }
ggrange <- function(output_wvi, Y) {
  if(!is.null(dim(Y))){
    if(length(dim(Y))>1){
      stop("Y must be a numeric vector, not a matrix or data frame.")
    }
  }
  particles <- output_wvi$particles
  K <- nrow(particles)
  ka <- c()
  cla <- c()
  meanxa <- c()
  minxa <- c()
  maxxa <- c()
  for (k in 1:K) {
    configk <- particles[k, ]
    nc <- max(configk)
    for (c in 1:nc) {
      ka <- c(ka, k)
      cla <- c(cla, c)
      meanxa <- c(meanxa, mean(Y[configk == c]))
      minxa <- c(minxa, min(Y[configk == c]))
      maxxa <- c(maxxa, max(Y[configk == c]))
    }
  }
  df <- data.frame(
    meanx = meanxa, minx = minxa, maxx = maxxa,
    cluster = as.factor(cla), k = ka
  )
  ggplot2::ggplot(df) +
    ggplot2::geom_point(ggplot2::aes(x = factor(.data$k), y = .data$meanx, col = .data$cluster)) +
    ggplot2::geom_errorbar(ggplot2::aes(x = factor(.data$k), ymin = .data$minx, ymax = .data$maxx, col = .data$cluster)) +
    ggplot2::labs(x = "particles", y = "mean and range of data within clusters")
}

#' Plot mean/range plot and data histogram side-by-side
#'
#' Shows ggrange plot and a histogram of the data together.
#'
#' @param output_wvi A list containing at least elements \code{particles} (matrix) from WASABI output.
#' @param Y Numeric vector of data.
#' @return If gridExtra is available, NULL (the plots are drawn). Otherwise, a list of plot objects.
#' @export
#' @examples
#' \dontrun{
#' ggrange_hist(output_wvi, Y)
#' }
ggrange_hist <- function(output_wvi, Y) {

  p1 <- ggrange(output_wvi, Y) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "left")

  p2 <- ggplot2::ggplot(data.frame(Y = Y)) +
    ggplot2::geom_histogram(ggplot2::aes(x = .data$Y)) +
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::theme_bw()

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(
      p1, p2,
      nrow = 1,
      layout_matrix = matrix(c(1, 1, 1, 2), nrow = 1),
      top = "Range and mean of data within clusters, with histogram of data"
    )
  } else {
    warning("gridExtra not installed; returning plots as a list.")
    return(list(plot1 = p1, plot2 = p2))
  }
}

#' Scatterplot: clusters vs 1d data, with jitter
#'
#' Makes a scatter plot of data vs cluster assignments, adding random vertical jitter for visual clarity.
#'
#' @param cls Numeric or factor vector of cluster assignments for each point.
#' @param Y Numeric vector of data.
#' @return A ggplot2 plot object.
#' @export
#' @examples
#' \dontrun{
#' ggscatter(output_wvi$particles[1,], Y)
#' }
ggscatter <- function(cls, Y) {
  if(!is.null(dim(Y))){
    if(length(dim(Y)) > 1){
      stop("Y must be a numeric vector, not a matrix or data frame.")
    }
  }
  n <- ifelse(is.null(dim(Y)), length(Y), nrow(Y)) 
  df_true <- data.frame(
    x = as.numeric(Y),
    cl = as.numeric(cls),
    jit = 0.2 * stats::rnorm(n)
  )
  ggplot2::ggplot(df_true) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$x,
        y = .data$cl + .data$jit,
        color = as.factor(.data$cl)
      )
    )
}

#' Grid of scatterplots: clusters in multiple particles (1D)
#'
#' Make a faceted scatter plot: each facet shows cluster assignments for one particle.
#'
#' @param output_wvi WASABI output; must contain \code{particles} and \code{part.weights}.
#' @param Y Numeric vector of data.
#' @return A ggplot2 plot object.
#' @export
#' @examples
#' \dontrun{
#' ggscatter_grid(output_wvi, Y)
#' }
ggscatter_grid <- function(output_wvi, Y) {
  if(!is.null(dim(Y))){
    if(length(dim(Y)) > 1){
      stop("Y must be a numeric vector, not a matrix or data frame.")
    }
  }
  n <- ifelse(is.null(dim(Y)), length(Y), nrow(Y)) 
  
  particles <- output_wvi$particles
  weights <- round(output_wvi$part.weights, 2)
  y <- 0.2 * stats::rnorm(n)
  df <- data.frame(
    x = NULL,
    cl = NULL,
    jit = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:length(weights)) {
    df_tmp <- data.frame(
      x = as.numeric(Y),
      cl = particles[i, ],
      jit = y,
      ID = i,
      w = weights[i]
    )
    df <- rbind(df, df_tmp)
  }

  ggplot2::ggplot(df) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$x,
        y = .data$cl + .data$jit,
        color = as.factor(.data$cl)
      )
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$ID))
}

#' Summary barplot for WASABI particles
#'
#' Plots both particle weights and corresponding number of clusters as a barplot for each sampled particle.
#'
#' @param output_wvi WASABI output; must contain \code{particles} and \code{part.weights}.
#' @param title Optional plot title.
#' @param legend Position of legend ("right" [default] or "bottom").
#' @return A ggplot2 barplot object.
#' @export
#' @examples
#' \dontrun{
#' ggsummary(output_wvi)
#' }
ggsummary <- function(output_wvi, title = NULL, legend = "right") {
  K <- length(output_wvi$part.weights)
  particles <- as.factor(seq_len(K))
  clusters <- apply(output_wvi$particles, 1, function(x) length(unique(x)))
  weights <- output_wvi$part.weights
  dat <- data.frame(
    particles = particles,
    Weight = weights,
    Number_of_clusters = clusters
  )

  # Manual pivot_longer using base R
  dat_long <- rbind(
    data.frame(
      particles = dat$particles,
      type = "Weight",
      value = dat$Weight,
      stringsAsFactors = FALSE
    ),
    data.frame(
      particles = dat$particles,
      type = "Number of clusters",
      value = dat$Number_of_clusters,
      stringsAsFactors = FALSE
    )
  )
  dat_long$type <- factor(dat_long$type, levels = c("Weight", "Number of clusters"))

  p <- ggplot2::ggplot(dat_long) +
    ggplot2::geom_bar(
      ggplot2::aes(
        x = .data$particles,
        y = .data$value,
        fill = .data$particles
      ),
      stat = "identity", width = 0.7
    ) +
    ggplot2::labs(x = "particles", y = "") +
    ggplot2::xlab("") +
    ggplot2::scale_fill_brewer(palette = "Accent") +
    ggplot2::facet_wrap(ggplot2::vars(.data$type), scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  
  if (legend == "bottom") {
    p <- p + ggplot2::theme(legend.position = "bottom")
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  return(p)
}


# ggsummary <- function(output_wvi, title = NULL, legend = "right") {
#   if (!requireNamespace("tidyr", quietly = TRUE)) {
#     stop("The 'tidyr' package is required for this function. Please install it.")
#   }
#   ind <- sort(output_wvi$part.weights, index.return = TRUE, decreasing = TRUE)
#   K <- length(output_wvi$part.weights)
#   particles <- as.factor(seq_len(K))
#   clusters <- apply(output_wvi$particles[ind$ix, ], 1, function(x) length(unique(x)))
#   weights <- output_wvi$part.weights[ind$ix]
#   dat <- data.frame(
#     particles = particles, Weight = weights,
#     Number_of_clusters = clusters
#   )
#   dat <- tidyr::pivot_longer(
#     dat,
#     cols = c("Weight", "Number_of_clusters"),
#     names_to = "type", values_to = "value"
#   )
#   dat$type <- ifelse(dat$type == "Number_of_clusters", "Number of clusters", dat$type)
#   dat$type <- factor(dat$type, levels = c("Weight", "Number of clusters"))

#   p <- ggplot2::ggplot(dat) +
#     ggplot2::geom_bar(
#       ggplot2::aes(
#         x = .data$particles,
#         y = .data$value,
#         fill = .data$particles
#       ),
#       stat = "identity", width = 0.7
#     ) +
#     ggplot2::labs(x = "particles", y = "") +
#     ggplot2::xlab("") +
#     ggplot2::scale_fill_brewer(palette = "Accent") +
#     ggplot2::facet_wrap(ggplot2::vars(.data$type), scales = "free") +
#     ggplot2::theme_bw() +
#     ggplot2::theme(axis.title.x = ggplot2::element_blank())
#   if (legend == "bottom") {
#     p <- p + ggplot2::theme(legend.position = "bottom")
#   }
#   if (!is.null(title)) {
#     return(p + ggplot2::ggtitle(title))
#   } else {
#     return(p)
#   }
# }

#' 2D cluster scatterplot
#'
#' Plots a scatter plot of the first two columns of feature data, colored by clusters.
#'
#' @param cls Numeric or factor vector of cluster assignments for each row.
#' @param Y Numeric matrix or data.frame with at least two columns.
#' @return A ggplot2 plot object.
#' @export
#' @examples
#' \dontrun{
#' ggscatter2d(output_wvi$particles[1,], Y)
#' }
ggscatter2d <- function(cls, Y) {
  if(is.null(dim(Y))){
    stop("Data must be 2-dimensional.")
  } else {
    if(length(dim(Y))==1){
      stop("Data must be 2-dimensional.")
    } else {
      if(length(dim(Y))>2){
        warning("Only first two dimensions used for plotting. Consider using ggscatter_pca instead.")
      }
    }
  }
  df_true <- data.frame(
    x1 = Y[, 1],
    x2 = Y[, 2],
    Cluster = as.factor(cls)
  )
  ggplot2::ggplot(df_true) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$x1,
        y = .data$x2,
        color = .data$Cluster
      )
    ) +
    ggplot2::theme_bw()
}

#' Grid of 2D cluster scatterplots from multiple particles
#'
#' For each particle, plot its clustering in the PC1-PC2 space or first-two-feature space, in a facet.
#'
#' @param output_wvi WASABI output; must contain \code{particles} and \code{part.weights}.
#' @param Y Numeric matrix/data.frame with at least two columns.
#' @param title Optional title for the plot.
#' @param legend Where to put the legend ("right" or "bottom").
#' @return A ggplot2 multi-facet plot object.
#' @export
#' @examples
#' \dontrun{
#' ggscatter_grid2d(output_wvi, Y)
#' }
ggscatter_grid2d <- function(output_wvi, Y, title = NULL, legend = "right") {
  if(is.null(dim(Y))){
    stop("Data must be 2-dimensional.")
  } else {
    if(length(dim(Y))==1){
      stop("Data must be 2-dimensional.")
    } else {
      if(length(dim(Y))>2){
        warning("Only first two dimensions used for plotting. Consider using ggscatter_pca instead.")
      }
    }
  }
  particles <- output_wvi$particles
  weights <- round(output_wvi$part.weights, 2)
  df <- data.frame(
    x1 = NULL,
    x2 = NULL,
    Cluster = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:nrow(particles)) {
    df_tmp <- data.frame(
      x1 = Y[, 1],
      x2 = Y[, 2],
      Cluster = particles[i, ],
      ID = paste0("Particle ", i, ", w = ", weights[i]),
      w = weights[i]
    )
    df <- rbind(df, df_tmp)
  }
  df$Cluster <- as.factor(df$Cluster)

  if(length(unique(df$Cluster))<=25){
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = .data$x1,
          y = .data$x2,
          color = .data$Cluster,
          shape = .data$Cluster
        )
    ) + scale_shape_manual(values = c(1:25))
  } else {
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = .data$x1,
          y = .data$x2,
          color = .data$Cluster
        )
      )
  }

  p <- p +
    ggplot2::ylab(expression("x"[2])) +
    ggplot2::xlab(expression("x"[1])) +
    ggplot2::facet_wrap(ggplot2::vars(.data$ID)) +
    ggplot2::theme_bw()
  if (legend == "bottom") {
    p <- p + ggplot2::theme(legend.position = "bottom")
  }
  if (!is.null(title) && title != FALSE) {
    p <- p + ggplot2::ggtitle(title)
  }
  return(p)
}

#' Grid of 2D cluster scatterplots from PCA coordinates
#'
#' Like \code{ggscatter_grid2d}, but uses principal component axes for the x/y axes.
#'
#' @param output_wvi WASABI output; must contain \code{particles} and \code{part.weights}.
#' @param Y Numeric matrix of dimension $n$ by $d$, with $d > 2$, containing the data points.
#' @param title Optional plot title.
#' @return A ggplot2 multi-facet plot object.
#' @export
#' @examples
#' \dontrun{
#' ggscatter_pca(output_wvi, Y)
#' }
ggscatter_pca <- function(output_wvi, Y, title = NULL) {
  if(is.null(dim(Y))){
    stop("Data must be d-dimensional with d > 2.")
  } else {
    if(length(dim(Y))==1){
      stop("Data must be d-dimensional with d > 2.")
    } 
  }

  out_pc <- stats::prcomp(Y)
  Ypc <- out_pc$x
  particles <- output_wvi$particles
  weights <- round(output_wvi$part.weights, 2)
  df <- data.frame(
    x1 = NULL,
    x2 = NULL,
    Cluster = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:nrow(particles)) {
    df_tmp <- data.frame(
      x1      = Ypc[, 1],
      x2      = Ypc[, 2],
      Cluster = particles[i, ],
      ID      = paste0("Particle ", i, ", w = ", weights[i]),
      w       = weights[i]
    )
    df <- rbind(df, df_tmp)
  }
  df$Cluster <- as.factor(df$Cluster)

  if(length(unique(df$Cluster))<=25){
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = .data$x1,
          y = .data$x2,
          color = .data$Cluster,
          shape = .data$Cluster
        )
    ) + scale_shape_manual(values = c(1:25))
  } else {
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = .data$x1,
          y = .data$x2,
          color = .data$Cluster
        )
      )
  }

  p <- p + 
    ggplot2::ylab(expression("PC"[2])) +
    ggplot2::xlab(expression("PC"[1])) +
    ggplot2::facet_wrap(ggplot2::vars(.data$ID)) +
    ggplot2::theme_bw()

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  return(p)
}

