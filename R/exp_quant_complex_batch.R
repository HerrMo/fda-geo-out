### batchtools experiment setup

library(batchtools)

packs <- c("dbscan", "fdaoutlier", "pROC", "data.table")

unlink("batch_exp/exp_fda_out", recursive = TRUE)
regis <- makeExperimentRegistry("batch_exp/exp_fda_out",
                                packages = packs,
                                source = c("R/data_generating_processes.R",
                                           "R/methods_and_performance_wrappers.R"))



# 1 Synthetic problems -------------------------------------------------------------------------------------------------

meths <- c("w_lof", "mds_lof", "do", "tv", "md")
mods <- paste("cDGP", c("fab1", "fab2", "fdao1", "fdao2"), sep = "_")

data_wrapper <- function(data, job, ...) {
  out_data(...)
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

res <- reduceResultsDataTable(fun = round, digits = 2)

pars <- unwrap(getJobPars())
res <- ijoin(pars, res)
head(res)

# save(res, file = paste0("data/res_complex_batch_", Sys.Date(), ".RData"))

