### batchtools experiment setup

library(batchtools)

# Experiment setup

packs <- c("dbscan", "pROC", "data.table", "classiFunc", "manifun")

unlink("batch_exp/exp_arrow", recursive = TRUE)
regis <- makeExperimentRegistry("batch_exp/exp_arrow",
                                packages = packs,
                                source = c("R/data_generating_processes.R",
                                           "R/methods_and_performance_wrappers.R"))



# Data -------------------------------------------------------------------------------------------------

# the output of data wrapper consists of
#   fun_obs: matrix with functional observations in rows
#   fun_dists: dist matrix of functional observations
#   out_labels: observation labels (0: inclass, 1: outlier)
#   out_indices: indices of outliers

arr_te <- read.table("data/raw/ArrowHead_TEST.txt")
arr_tr <- read.table("data/raw/ArrowHead_TRAIN.txt")

arr <- rbind(arr_te, arr_tr)

data_wrapper <- function(data, job, dis = dis, r_out = r_out, ...) {

  dat_in <- data[data$V1 == 0, ]
  dat_out <- data[data$V1 != 0, ]

  # remove label independent structural outliers
  # - in all three classes there a observations which can be considered structural anomalies
  # - need to be removed if performance assessment is based on class label information
  dat_in <- dat_in[-c(17, 40, 45), ]
  dat_outs <- dat_out[-c(65, 81, 93, 122), ]

  n_in <- nrow(dat_in)
  n_out <- nrow(dat_out)

  dat_in <- as.matrix(dat_in[, -1])
  dat_out <- as.matrix(dat_out[, -1])

  n_out_smpl <- round(n_in * r_out)
  i_out_smpl <- sample(seq_len(n_out), n_out_smpl)

  dat_out_smpl <- dat_out[i_out_smpl, ]

  data <- list(data = rbind(dat_in, dat_out_smpl),
               true_outliers = seq_len(n_out_smpl) + n_in)

  if (dis == "l2") {
    fun_dists <- dist(data$data)
  } else if (dis == "dtw") {
    fun_dists <- as.matrix(classiFunc::computeDistMat(data$data, method = "rucrdtw"))
  } else {
    grid <- seq(0, 1, length.out = ncol(data$data))
    fun_dists <- as.matrix(wasser_dist(data$data, grid = grid))
  }
  out_labels <- seq_len(nrow(data$data)) %in% data$true_outliers

  dat <- list(fun_obs = data$data,
              fun_dists = fun_dists,
              out_labels = out_labels + 0,
              out_indices = data$true_outliers)
  dat
}

addProblem("arrowhead", data = arr, fun = data_wrapper, seed = 2)

pdes = list(
  arrowhead =
    CJ(
      r_out = c(0.1, 0.05),
      dis = c("l2", "dtw", "wasser")
    )
)


# 2 Methods and performance --------------------------------------------------------------------------------------------
meths <- c("mds_lof")

performance <- function(dat, scores) {
  ranks <- rank(-scores)
  out_ind <- dat$out_indices
  out_lbls <- dat$out_labels

  data.table(
    auc = auc(roc(out_lbls, scores, direction = "<")),
    tpr = tpr(ranks, out_ind))
}

mds_lof_wrapper <- function(data, job, instance, meth, ...) {
  scores <- mds_lof(instance)

  performance(instance, scores)
}

# not used
w_lof_wrapper <- function(data, job, instance, meth, ...) {
  scores <- w_lof(instance)

  performance(instance, scores)
}

d_lof_wrapper <- function(data, job, instance, meth, ...) {
  d_dist <- function(mat, grid) {

    derivs <- t(apply(mat, 1, function(f) my_deriv(grid = grid, vals = f)))

    dfuns <- dist(mat)
    dders <- dist(derivs)

    as.matrix(dders)
  }

  my_deriv <- function(grid, vals){
    diff(vals)/diff(grid)
  }

  grid <- seq(0, 1, length.out = ncol(instance$fun_obs))
  ddists <- d_dist(instance$fun_obs, grid = grid)

  emb <- cmdscale(ddists, k = 5)
  scores <- lof(emb, minPts = 0.75 * nrow(emb))

  performance(instance, scores)
}

do_wrapper <- function(data, job, instance, meth, ...) {
  scores <- do(instance)

  performance(instance, scores)
}

tv_wrapper <- function(data, job, instance, meth, ...) {
  scores <- tv(instance)

  performance(instance, scores)
}

md_wrapper <- function(data, job, instance, meth, ...) {
  scores <- md(instance)

  scrs <- apply(scores, 2, performance, dat = instance)
  m_scrs <- do.call(rbind, scrs)
  rownames(m_scrs) <- names(scrs)
  m_scrs
}


for (met in meths) {
  meth <- paste0(met, "_wrapper")
  addAlgorithm(paste0(met, "_wrapper"), fun = get(meth))
}

addExperiments(pdes, repls = 500)
summarizeExperiments()

testJob(1)


regis$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 6)

submitJobs()

getStatus()


res <- reduceResultsDataTable(fun = round, digits = 2)

pars <- unwrap(getJobPars())
res_arrow <- ijoin(pars, res)
head(res_arrow)

save(res_arrow, file = paste0("data/res_arrow_batch_", Sys.Date(), ".RData"))

