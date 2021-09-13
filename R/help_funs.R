# dist and lof rank calculation helpers

dists <- function(x) as.matrix(proxy::dist(x))
frank <- function(x, ndim = 2, k = 0.75 * nrow(extract_points(x))) {rank(-lof(extract_points(x, ndim = ndim), minPts = k))}


# simulate data from Hernandez and Munoz ------------------------------------------------------------------------------
hm_dat <- function(n_in = 95, n_out = 5) {
  a <- rnorm(n_in, 5, 2)
  b <- rnorm(n_out, 5, sqrt(3))

  t <- 1:500
  xt <- t/500

  inclass <- t(simplify2array(lapply(a, function(x) {x + 0.05 * t + sin(pi * (xt)^2)})))
  outclass <- t(simplify2array(lapply(b, function(x) {x + 0.05 * t + cos(20 * pi * xt)})))
  rbind(inclass, outclass)
}

# wasserstein dist and help funs --------------------------------------------------------------------------------------
#' @export
wasserstein_metric <- function(f1, f2, grid, a = 1, norm = FALSE) {
  # f1 = mu0, f2 = mu1 in the paper

  require(caTools) # for trapz (numerical integration)
  require(transport) # for wasserstein1d

  if (norm) return(wasserstein1d(f1, f2)) # computes standard wasserstein-1 distance

  # compute standard measure of domain
  omega <- tail(grid, 1) - grid[1]

  f_diff <- f2 - f1
  c <- trapz(grid, f_diff) * (1/omega)

  mu <- f_diff - c
  if (sum(mu) == 0) return(0) # saves some time if f1 = f2

  mu1 <- (mu + abs(mu)) / 2 # compute positive part
  mu2 <- (mu - abs(mu)) / 2 # compute negative part

  wasserstein1d(mu1, mu2) + omega/a * abs(c)
}

#' @export
wasser_dist <- function(dat, ...) {
  nfun <- nrow(dat)
  combs <- combn(1:nfun, 2)
  m_dist <- diag(x = 0, nfun)
  for (i in seq_len(ncol(combs))) {
    ind <- combs[, i]
    ind1 <- ind[1]
    ind2 <- ind[2]
    m_dist[ind1, ind2] <- wasserstein_metric(dat[ind1, ],
                                             dat[ind2, ],
                                             ...)
  }
  m_dist[lower.tri(m_dist)] <- t(m_dist)[lower.tri(m_dist)]
  m_dist
}

d_dist <- function(mat, a = 0.5, grid, d_only = TRUE) {

  derivs <- t(apply(mat, 1, function(f) my_deriv(grid = grid, vals = f)))

  dfuns <- dist(mat)
  dders <- dist(derivs)

  a * as.matrix(dfuns) + (1 - a) * as.matrix(dders)
}

my_deriv = function(grid, vals){
  diff(vals)/diff(grid)
}


# help fun S3 class to extract embedding coordinates ------------------------------------------------------------------
#' @export
extract_points <- function(x, ...) {
  UseMethod("extract_points")
}

# S3 method for isomap
#' @export
extract_points.isomap <- function(embedding, ndim = dim(embedding$points)[2]) {
  embedding$points[, 1:ndim, drop = FALSE]
}

# S3 method for umap
#' @export
extract_points.umap <- function(embedding, ndim = dim(embedding$layout)[2]) {
  embedding$layout[, 1:ndim, drop = FALSE]
}

# S3 method for matrix output, e.g. mds
#' @export
extract_points.matrix <- function(embedding, ndim = dim(embedding)[2]) {
  embedding[, 1:ndim, drop = FALSE]
}


# embedding wrapper ---------------------------------------------------------------------------------------------------
embed <- function(dist_mat, method = c("isomap", "umap", "mds"), ...) {
  method <- match.arg(method, c("isomap", "umap", "mds"))

  if (inherits(dist_mat, "dist")) dist_mat <- as.matrix(dist_mat)

  emb <-
    switch(
      method,
      "isomap" = vegan::isomap(dist_mat, ...),
      "umap" = umap::umap(dist_mat, input = "dist", ...),
      "mds" = cmdscale(dist_mat, ...)
    )

 emb
}



