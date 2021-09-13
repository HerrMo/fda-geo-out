### batchtools experiment setup

library(batchtools)

packs <- c("dbscan", "fdaoutlier", "pROC", "data.table")

# unlink("batch_exp/exp_l2vl10", recursive = TRUE)
regis <- makeExperimentRegistry("batch_exp/exp_l2vl10",
                                packages = packs,
                                source = c("R/data_generating_processes.R",
                                           "R/methods_and_performance_wrappers.R"))



# 1 Synthetic problems -------------------------------------------------------------------------------------------------

meths <- c("mds_lof")
mods <- 2

data_wrapper <- function(data, job, mod = mod, dis = dis, n = n, r_out = r_out, ...) {
  dat <- get(paste0("simulation_model", mod))(n = n, p = 50, outlier_rate = r_out, plot = FALSE)

  if (dis == "l2") {
    fun_dists <- dist(dat$data)
  } else {
    fun_dists <- dist(dat$data, method = "minkowski", p = 10)
  }
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
      dis = c("l2", "l10"),
      r_out = c(0.1, 0.05),
      mod = mods
    ),
  n1000 =
    CJ(
      n = 1000,
      dis = c("l2", "l10"),
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


res <- reduceResultsDataTable(fun = round, digits = 2)

pars <- unwrap(getJobPars())
res_l2vsl10 <- ijoin(pars, res)
head(res_l2vsl10)

save(res_l2vsl10, file = paste0("data/res_l2vsl10_batch_", Sys.Date(), ".RData"))

