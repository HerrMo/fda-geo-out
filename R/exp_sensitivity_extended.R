### Sensitivity analysis

# source(here::here("R/setup.R")) # uncomment setup.R for stand-alone usage

dim_cor_frank <- function(points, dims, minPts = 0.75 * nrow(points)) {
  franks <- vapply(dims, function(dim) frank(points, ndim = dim, k = minPts), FUN.VALUE = numeric(nrow(points)))
  cor_mat <- cor(franks, method = "spearman")
  up_tri <- upper.tri(cor_mat)
  cor_vec <- cor_mat[up_tri]
  dim_grid <- expand.grid(dims, dims)
  entry_nams <- matrix(paste(dim_grid[, 1], dim_grid[, 2], sep = "vs"), nrow = length(dims))
  names(cor_vec) <- entry_nams[up_tri]
  cor_vec
}


## settings

dat_sets <-  c("ECG", "Octane", "Weather", "Tecator", "Wine")
n_datsets <- length(dat_sets)
dims <- 1:20
max_emb_dim <- max(dims)

metrics <- c("l.5", "l1", "l2", "l4", "l10", "wasser")

for (dat in dat_sets) {
  load(here::here(paste0("data/processed/", dat, ".RData")))
}

dat_names <- c("out_ecg200", "oct", "tec", "wine", "weather")

dist_lists <- vector("list", length(dat_names))
names(dist_lists) <- dat_names
for (dat_set in dat_names) {
  dat <- get(dat_set)
  dat <- as.matrix(dat)
  dist_lists[[dat_set]] <-
    lapply(metrics, function(metric) {
      dis <- switch(metric,
                    "l.5" = dist(dat, "minkowski", p = 0.5),
                    "l1" = dist(dat, "minkowski", p = 1),
                    "l2" = dist(dat),
                    "l4" = dist(dat, "minkowski", p = 4),
                    "l10" = dist(dat, "minkowski", p = 10),
                    "wasser" = wasser_dist(dat,
                                         grid = seq(0, 1, length.out = ncol(dat))))
      dis
    })
  names(dist_lists[[dat_set]]) <- metrics
}


# embedding dimensionality: GOF ----------------------------------------------------------------------------------------
names(dist_lists) <- dat_sets

l_expl_vari <- l_embs_k20 <- vector("list", length(metrics))
names(l_expl_vari) <- names(l_embs_k20) <-  metrics
for (met in metrics) {
  # list of distances by metric
  dist_list <- lapply(dist_lists, function(dat) dat[[met]])
  # mds embeddings for with dim = max_emb_dim
  mds_embs_k20 <- lapply(dist_list, function(d) cmdscale(d, k = max_emb_dim, eig = TRUE))
  l_embs_k20[[met]] <- mds_embs_k20
  par_grid <- expand.grid(data = dat_sets, dims = dims)

  expl_vari <- mapply(
    function(dat, dim) {
      evals <- mds_embs_k20[[dat]]$eig
      pos_evals <- evals[evals > 0]
      sum(pos_evals[1:dim]) / sum(pos_evals) # Mead version
      # sum(abs(evals[1:dim])) / sum(abs(evals)) # Mardia version
    },
    dat = par_grid$data,
    dim = par_grid$dims)

  l_expl_vari[[met]] <- data.table(
    "expl_variation" = expl_vari,
    par_grid
  )
}

# embedding dimensionality: LOF rank correlation -----------------------------------------------------------------------
cor_dims <- 2:20
n_cor_dims <- length(cor_dims)

# number of correlations (without autocorrelations), i.e. number of cells in upper/lower triangle of the correlation matrix
out_len <- (n_cor_dims * (n_cor_dims - 1)) / 2

l_cor_mats <- vector("list", length(metrics))
names(l_cor_mats) <- metrics

for (met in metrics) {
  l_cor_mats[[met]] <-
    vapply(l_embs_k20[[met]],
           function(x) dim_cor_frank(x$points, dims = cor_dims), FUN.VALUE = numeric(out_len))
}

l_cor_dt <- lapply(l_cor_mats, as.data.table, keep.rownames = TRUE)
dt_cors <- rbindlist(l_cor_dt, idcol = "file")
setnames(dt_cors, c("file", "rn"), c("Metric", "Dims"))
num_cols <- dat_sets
dt_cors[, (num_cols) := round(.SD, 2), .SDcols = num_cols]

cmb_of_interest <- c("2vs5", "5vs20")
dt_cors <- dt_cors[Dims %in% cmb_of_interest, ]

