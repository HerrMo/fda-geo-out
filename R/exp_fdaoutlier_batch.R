### batchtools experiment setup

library(batchtools)

# Experiment setup

packs <- c("dbscan", "fdaoutlier", "pROC", "data.table")

unlink("batch_exp/exp_fdaoutlier_deriv", recursive = TRUE)
regis <- makeExperimentRegistry("batch_exp/exp_fdaoutlier_deriv",
                                packages = packs,
                                source = c("R/data_generating_processes.R",
                                           "R/methods_and_performance_wrappers.R"))



# 1 Synthetic problems -------------------------------------------------------------------------------------------------

meths <- c("w_lof", "mds_lof", "do", "tv", "md")
mods <- c(1:9)

data_wrapper <- function(data, job, mod = mod, n = n, r_out = r_out, ...) {
  dat <- get(paste0("simulation_model", mod))(n = n, p = 50, outlier_rate = r_out, plot = FALSE)

  fun_dists <- dist(dat$data)
  out_labels <- seq_len(n) %in% dat$true_outliers

  dat <- list(fun_obs = dat$data,
              fun_dists = fun_dists,
              out_labels = out_labels + 0,
              out_indices = dat$true_outliers)
  dat
}


for (i in c(100, 1000)) {
  addProblem(paste0("n", i), fun = data_wrapper)
}

pdes = list(
  n100 =
    CJ(
      n = 100,
      r_out = c(0.1, 0.05),
      mod = mods
    ),
  n1000 =
    CJ(
      n = 1000,
      r_out = 0.01,
      mod = mods
    )
)


# 2 Algos --------------------------------------------------------------------------------------------------------------
performance <- function(dat, scores) {
  ranks <- rank(-scores)
  out_ind <- dat$out_indices
  out_lbls <- dat$out_labels

  data.table(
    auc = auc(roc(out_lbls, scores, direction = "<")),
    tpr = tpr(ranks, out_ind))
}

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

mds_lof_wrapper <- function(data, job, instance, meth, ...) {
  scores <- mds_lof(instance)

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

# gather results
# res = do.call(rbind, reduceResultsList(fun = round, digits = 2))
# res$job.id <- seq_len(nrow(res))
# head(res)

res_deriv <- reduceResultsDataTable(fun = round, digits = 2)

pars <- unwrap(getJobPars())
res_deriv <- ijoin(pars, res_deriv)
head(res_fdaoutlier)

save(res_fdaoutlier, file = paste0("data/res_fdaoutlier_batch_", Sys.Date(), ".RData"))

