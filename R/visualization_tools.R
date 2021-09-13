#' @export
plot_funs <- function(data, ...) {
  UseMethod("plot_funs")
}

#' @export
# default function for data in fundata in matrix format
plot_funs.default <- function(data, col = NULL, args = NULL) {
  n <- nrow(data)
  grid_len <- ncol(data)
  df_dat <- data.frame(
    args = if (is.null(args)) rep(1:grid_len, n) else args,
    vals = c(t(data)),
    id = as.factor(rep(1:n, each = grid_len))
  )

  if (!is.null(col)) df_dat$col <- rep(col, each = grid_len)

  ggplot(df_dat) +
    geom_line(aes(x = args,
                  y = vals,
                  group = id,
                  colour = if (is.null(col)) {id} else {col})) +
    theme(legend.position = "None")
}

# function to plot embeddings
#' @export
plot_emb <- function(embedding, ...) {
  UseMethod("plot_emb")
}

# default method for embedding data in 2d matrix format
#' @export
plot_emb.default <- function(pts, color = NULL, size = 1, ...) {

  dat <- data.frame(dim1 = pts[, 1],
                    dim2 = pts[, 2],
                    color = 1:nrow(pts))

  if (!is.null(color)) dat$color <- color

  p <- ggplot(dat) +
    geom_point(aes(x = dim1,
                   y = dim2,
                   colour = color),
               size = size) +
    theme(legend.position = "Non") +
    ggtitle(label = "2d-embedding")
  p
}

# for embeddings coordinates in matrix format
#' @export
plot_emb.matrix <- function(embedding, color = NULL, labels_off = TRUE, labels = NULL, size = 1, ...) {
  # TODO argument checking (min 2-d data, etc)

  pts <- extract_points(embedding, 2)
  p <- plot_emb.default(pts, color = color, labels = labels, size = size, ...)
  if (!labels_off) p <- if (is.null(labels)) {
    p + ggrepel::geom_text_repel(aes(x = dim1, y = dim2, label = 1:nrow(pts)), size = label_size)
  } else {
    p + ggrepel::geom_text_repel(aes(x = dim1, y = dim2, label = labels), size = label_size)
  }
  p
}

plot_emb_temp <- function(pts, color = NULL, size = 1, pch = rep(1, nrow(pts)), ...) {

  dat <- data.frame(dim1 = pts[, 1],
                    dim2 = pts[, 2],
                    color = 1:nrow(pts),
                    pch = as.factor(pch))

  if (!is.null(color)) dat$color <- color

  p <- ggplot(dat) +
    geom_point(aes(x = dim1,
                   y = dim2,
                   colour = color,
                   pch = pch),
               size = size) +
    theme(legend.position = "Non") +
    ggtitle(label = "2d-embedding")
  p
}
