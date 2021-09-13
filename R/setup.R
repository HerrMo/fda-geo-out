## setup for standalone usage of the experiment files exp_real_qual.R and exp_sensitivity_extended.R

## packages

if (Sys.info()["user"] != "fabians") {
  devtools::load_all("~/projects/manifun")
  devtools::load_all("~/projects/simfun")
} else {
  devtools::load_all("~/fda/manifun")
  devtools::load_all("~/fda/simfun")
  plot_funs.matrix <- manifun::plot_funs.default
  plot_funs.data.frame <- manifun::plot_funs.default
  plot_emb.matrix <- manifun::plot_emb.matrix
}

library(data.table)
library(patchwork)
library(fda)
library(fda.usc)
library(parallel)
library(dbscan)
library(ggplot2)
library(igraph)
library(StatPerMeCo)
library(pROC)
library(cowplot)
library(latex2exp)
library(fdaoutlier)
library(viridis)


## help functions

# dist calculation and lof ranks
dists <- function(x) as.matrix(proxy::dist(x))
frank <- function(x, ndim = 2, k = 0.75 * nrow(extract_points(x))) {rank(-lof(extract_points(x, ndim = ndim), minPts = k))}

## plotting funs
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

