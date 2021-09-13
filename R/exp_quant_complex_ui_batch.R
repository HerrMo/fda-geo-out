### batchtools experiment setup

library(batchtools)

packs <- c("dbscan", "fdaoutlier", "pROC", "data.table", "umap", "vegan")

unlink("batch_exp/exp_uimap", recursive = TRUE)
regis <- makeExperimentRegistry("batch_exp/exp_uimap",
                                packages = packs,
                                source = c("R/data_generating_processes.R",
                                           "R/methods_and_performance_wrappers.R"))



# 1 Synthetic problems -------------------------------------------------------------------------------------------------

meths <- c("imap_lof", "umap_lof")
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

imap_lof_wrapper <- function(data, job, instance, meth, ...) {
  emb <- vegan::isomap(as.matrix(instance$fun_dists), ndim = 5, ...)
  scores <- lof(emb$points, minPts = 0.75 * nrow(instance$fun_obs))

  performance(instance, scores)
}

umap_lof_wrapper <- function(data, job, instance, n_neighbors = n_neighbors, ...) {
  pars <- umap.defaults
  pars$n_neighbors <- n_neighbors
  pars$n_components <- 5
  pars$input <- "dist"
  emb <- umap::umap(as.matrix(instance$fun_dists), config = pars)
  scores <- lof(emb$layout, minPts = 0.75 * nrow(instance$fun_obs))

  performance(instance, scores)
}

ades <-  list(
  # imap_lof_wrapper =
  #   CJ(k = c(5, 90)
  #   ),
  umap_lof_wrapper =
    CJ(
      n_neighbors = c(5, 90)
    )
)

for (met in meths) {
  meth <- paste0(met, "_wrapper")
  addAlgorithm(paste0(met, "_wrapper"), fun = get(meth))
}

addExperiments(pdes, ades, repls = 500)
summarizeExperiments()

testJob(8001)


regis$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 6)

submitJobs()

getStatus()

res <- reduceResultsDataTable(fun = round, digits = 2)

pars <- unwrap(getJobPars())
res_umap <- ijoin(pars, res)
head(res_umap)

save(res_umap, file = paste0("data/res_umap_batch_", Sys.Date(), ".RData"))

