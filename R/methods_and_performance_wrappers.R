## performance

# wrapper for pROC::auc
# scores: outlier scores, the higher the more outlying
# out_labels: observation labels with 1 outlier and 0 normal observation

# TRUE positive outlier rate
# scores: outlier scores of all observations, the higher the more outlying
# true_outliers: the indices of the true outliers in the data set
tpr <- function(out_ranks, true_outliers) {
  n_out <- length(true_outliers)

  tp <- which(out_ranks <= n_out) %in% true_outliers
  tpr <- sum(tp)/n_out
  tpr
}


## methods
# wrapper for mds + lof
w_lof <- function(out_dat, minPts = 0.75 * nrow(out_dat$fun_dists)) {
  require(dbscan)

  scores <- lof(out_dat$fun_dists, minPts = minPts)
  scores
}

mds_lof <- function(out_dat, dim = 5, minPts = 0.75 * nrow(out_dat$fun_dists)) {
  require(dbscan)

  emb <- cmdscale(as.matrix(out_dat$fun_dists), k = dim)
  scores <- lof(emb, minPts = minPts)
  scores
}

# wrapper for fdaoutlier::dir_out
do <- function(out_dat) {
  require(fdaoutlier)

  scores <- dir_out(out_dat$fun_obs)$var_outlyingness
  scores
}

# wrapper for fdaoutlier::total_variational_depth
tv <- function(out_dat) {
  require(dbscan)

  # negative scores so that the rank orientations is the same for all scoring methods
  scores <- -total_variation_depth(out_dat$fun_obs)$tvd
  scores
}

# wrapper for fdaoutlier::muod
md <- function(out_dat) {
   require(fdaoutlier)

   scores <- muod(out_dat$fun_obs)$indices
   scores
}




