out_data <- function(mod, n, r_out, p = 50) {
  dat <- simulate_complex_data(dgp = mod, r_out = r_out, p = p, n = n)

  fun_dists <- dist(dat$data)
  out_labels <- seq_len(nrow(dat$data)) %in% dat$true_outliers

  dat <-
    list(fun_obs = dat$data,
         fun_dists = fun_dists,
         out_labels = out_labels + 0,
         out_indices = dat$true_outliers)
  dat
}


simulate_complex_data <- function(dgp, n, r_out, p) {
  # models:
  #  dgp1
  #  dgp2
  #  fdaoutlier simulation models

  require(fdaoutlier)

  n_out <- round(n * r_out)
  n_in <- n - n_out

  switch(dgp,
         cDGP_fab1 = dgp_fab(n_in = n_in, n_out = n_out, p = p, model = 1),
         cDGP_fab2 = dgp_fab(n_in = n_in, n_out = n_out, p = p, model = 2),
         cDGP_fdao1 = dgp_fdao1(n_in = n, r_out = r_out, p = p),
         cDGP_fdao2 = dgp_fdao2(n_in = n, r_out = r_out, p = p))
}


dgp_fdao1 <- function(n_in, r_out,  p) {
  require(fdaoutlier)

  mod1 <- simulation_model1(n = n_in, p = p, outlier_rate = r_out / 2, plot = FALSE)
  mod2 <- simulation_model8(n = n_in, p = p, outlier_rate = r_out / 2, plot = FALSE)
  mod2_out <- mod2$data[mod2$true_outliers, ]

  # randomly replace n_out normal class obs from mod1 with the outliers from mod2
  repl <- sample(which(!(seq_len(n_in) %in% mod1$true_outliers)), length(mod2$true_outliers))
  mod1$data[repl, ] <- mod2_out
  mod1$true_outliers <- sort(c(mod1$true_outliers, repl))
  mod1
}

dgp_fdao2 <- function(n_in, r_out, p) {
  require(fdaoutlier)

  mod1 <- simulation_model1(n = n_in, p = p, outlier_rate = 0, plot = FALSE)

  n_out <- round(n_in * r_out)

  out_mods <- 1:9
  out_mods <- sample(out_mods, n_out, replace = n_out > length(out_mods))

  out_mat <-
    t(vapply(out_mods,
             function(x) {
               temp_mod <- get(paste0("simulation_model", x))(n = 10, p = p, outlier_rate = 0.1, plot = FALSE)
               unname(temp_mod$data[temp_mod$true_outliers, ])
             },
             numeric(p)))

  mod1$data <- rbind(mod1$data, out_mat)
  mod1$true_outliers <- n_in + seq_len(n_out)
  mod1
}

dgp_fab <- function(n_in, n_out, p, model) {

  grid <- seq(0, 1, length.out = p)

  if (model == 1) {
    bs_df <- 15
    sd_noise <- .1
    inlier_warped <- replicate(n_in, pbeta(grid, runif(1, 4, 6), runif(1, 4, 6))) |> t()
    outlier_warped <- replicate(n_out, pbeta(grid, runif(1, 3, 4), runif(1, 3, 4)))  |> t()
  } else {
    bs_df <- 25
    sd_noise <- .15
    inlier_warped <- replicate(n_in, pbeta(grid, runif(1, 3, 8), runif(1, 3, 8))) |> t()
    outlier_warped <- 0.5 * replicate(n_out,
                                      pbeta(grid, runif(1, .1, 3), runif(1, .1, 3)) +
                                        pbeta(grid, runif(1, 3, 8), runif(1, 3, 8))) |> t()
  }

  bs_coefs <- scale(rnorm(bs_df))
  template <- function(x) splines::bs(x, df = bs_df) %*%  bs_coefs

  warped <- rbind(inlier_warped, outlier_warped)
  data <- apply(warped, 1, template) |> t()
  data <- data +
    matrix(rnorm(length(grid) * (n_in + n_out)), ncol = length(grid)) * sd_noise

  list(data = data, true_outliers = n_in + cumsum(rep(1, n_out)))
}

